#include <mpich/mpi.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/sysinfo.h> // number of cpu
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <cstdlib>
#include <functional>
#include <cassert>
#include <atomic>
//#include "timing.h"
#include <chrono>

#include "renderer.h"
#include "thread-pool.h"
#include "timing.h"

#define macroStringfy(x) #x

#define MSG_TAG_REGISTER 0
#define MSG_TAG_WORKER_RESULT 1
#define MSG_TAG_WORKER_JOB 2
#define MSG_TAG_WORKER_TILE 3

#define RECV_FREQ_MS 100

constexpr int MAX_ELEMENTS = 300;

#include "image.h"
#include "constants.h"

using std::cout;
using std::endl;


inline std::ostream& operator<<(std::ostream& os, const glm::vec3& v) {
    os << '[' << v.x << ',' << v.y << ',' << v.z << ']';
    return os;
}


void abort_msg(const char msg[])
{
    fprintf(stderr, "ABORT: %s\n", msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}



struct TileImage
{
    int tile_id;
    float data[TILE_SIZE*TILE_SIZE*3];

    static bool type_done;
    static MPI_Datatype createMPIType() {
        // TODO: This is hack.
        // TODO: It must be freed when not needed
        static MPI_Datatype TileImage_type;
        if (!type_done) {
            int lengths[2] = { 1, TILE_SIZE*TILE_SIZE*3 };
            const MPI_Aint displacements[2] = { 0, offsetof(TileImage, data) };
            MPI_Datatype types[2] = { MPI_INT, MPI_FLOAT };
            if (MPI_SUCCESS != MPI_Type_create_struct(2, lengths, displacements, types, &TileImage_type)) {
                abort_msg("MPI_Type_create_struct");
            }
            if (MPI_SUCCESS != MPI_Type_commit(&TileImage_type))            {
                abort_msg("MPI_Type_commit");
            }
            type_done = true;
        }
        return TileImage_type;
    }

    static void freeMPIType() {
        if (type_done) {
            MPI_Datatype TileImage_type = TileImage::createMPIType();
            MPI_Type_free(&TileImage_type);
            type_done = false;
        }
    }

};

bool TileImage::type_done = false;

void iterateImageTile(const TileInfo& tile, std::function<void(std::size_t cn, std::size_t rn)> callback) {
    for (std::size_t rn = tile.r_start; rn < tile.r_stop; rn++) {
        for (std::size_t cn = tile.c_start; cn < tile.c_stop; cn++) {
            // Sample position (0,0) is bottom left corner of image, but image.setColor(0,0) is
            // position of top left corner. We need to invert Y axis here.
            callback(cn, frame_rows - 1 - rn);
        }
    }
}


void test_feelImage(Image& image) {
    for (int i = 0; i < TileInfo::numberOfTiles(); i++) {
        iterateImageTile(i, [i, &image](std::size_t cn, std::size_t rn) {
            float f = (float)i / TileInfo::numberOfTiles();
            image.setColor(cn, rn, glm::vec3{f});
        });
    }
}


std::unique_ptr<TileImage> getImageRegion(const Image& image, int tile_id)
{
    TileInfo tile(tile_id);
    std::unique_ptr<TileImage> result = std::make_unique<TileImage>();
    result->tile_id = tile_id;
    std::size_t i = 0;
    iterateImageTile(tile, [&i, &image, result = result.get()](std::size_t cn, std::size_t rn) {
        auto c = image.getColor(cn, rn);
        result->data[3*i] = c.r;
        result->data[3*i + 1] = c.g;
        result->data[3*i + 2] = c.b;
        i++;
    });
    return result;
}


void test_iterateImageTile() {
    Image image{frame_columns, frame_rows};
    test_feelImage(image);
    image.writeToFileBMP("output.bmp");
}


void fillImageWithRegion(const TileImage& tile, Image& image) {
    std::size_t i = 0;
    iterateImageTile(tile.tile_id, [&](std::size_t cn, std::size_t rn) {
        glm::vec3 c{tile.data[3*i],
                    tile.data[3*i + 1],
                    tile.data[3*i + 2]};
        image.setColor(cn, rn, c);
        i++;
    });
}


void test_iterateImageTile2() {
    Image image{frame_columns, frame_rows};
    test_feelImage(image);
    Image image2{frame_columns, frame_rows};

    for (int i = 0; i < TileInfo::numberOfTiles()/2; i++) {
        auto r = getImageRegion(image, i);
        fillImageWithRegion(*r, image2);
    }
    image2.writeToFileBMP("output.bmp");
}


void workerSendImageTile(const Image& image, int tile_id) {
    auto tile_data = getImageRegion(image, tile_id);

    MPI_Datatype TileImage_type = TileImage::createMPIType();
    if (MPI_SUCCESS != MPI_Send(tile_data.get(), 1, TileImage_type, 0 /*master*/, MSG_TAG_WORKER_TILE, MPI_COMM_WORLD))
    {
        abort_msg("MPI_Send: workerSendResult");
    }
}



void receiveImageTile(Image& image) {
    constexpr int tag = MSG_TAG_WORKER_TILE;
    {
        int flags = 0;
        while(!flags) {
            if (MPI_SUCCESS != MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flags, MPI_STATUS_IGNORE)){
                abort_msg("receiveImageTile: MPI_Iprobe");
            }
            msleep(RECV_FREQ_MS);
        }
    }

    MPI_Status status;
    std::unique_ptr<TileImage> result = std::make_unique<TileImage>();
    MPI_Datatype TileImage_type = TileImage::createMPIType();
    if (MPI_SUCCESS != MPI_Recv(result.get(), 1, TileImage_type, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status)) {
        abort_msg("receiveImageTile: MPI_Recv");
    }

    cout << "Received tile: " << result->tile_id << endl;
    // Fill image by received data
    fillImageWithRegion(*result, image);
}




struct ImageGather {
    Image image{frame_columns, frame_rows};
    std::thread gather_thread;

    ImageGather() {
    }

    void gatherTiles() {
        gather_thread = std::thread{[this](){
            for (int i = 0; i < (int)TileInfo::numberOfTiles(); i++) {
                receiveImageTile(image);
            }
        }};
    }

    void join() {
        gather_thread.join();
    }
};




struct WorkerInfo
{
    //Worker(int worker_rank, int cpus) : worker_rank{worker_rank}, cpus{cpus} {}
    // Worker& operator=(const Worker& w) {
    //     wor
    // }

    std::string info()
    {
        std::stringstream ss;
        ss << "rank: " << worker_rank << " cpus: " << cpus;
        return ss.str();
    }
    // TODO: WOW: why I can not use const here? (providing that I have constructor)
    // test-distr.cc:27:7: note: ‘Worker& Worker::operator=(const Worker&)’ is implicitly deleted because the default definition would be ill-formed:
    // 27 | class Worker
    int worker_rank;
    int cpus; 
};

struct WorkerResult
{
    int worker_rank;
    int input;
};


void masterSendWorkerStop(int rank) {
    const int data = -1;
    printf("Mater: send job rank: %d STOP\n", rank);
    if (MPI_SUCCESS != MPI_Send(&data, 1, MPI_INT, rank, MSG_TAG_WORKER_JOB, MPI_COMM_WORLD))
    {
        abort_msg("MPI_Send");
    }
}

void masterSendWorkerJob(int rank, const std::vector<int> data) {
    //printf("Mater: send job rank: %d x=%d\n", const, data);
    if (MPI_SUCCESS != MPI_Send(data.data(), data.size(), MPI_INT, rank, MSG_TAG_WORKER_JOB, MPI_COMM_WORLD))
    {
        abort_msg("MPI_Send");
    }
}


void workerSendWorkDone(int tile_id) {
    if (MPI_SUCCESS != MPI_Send(&tile_id, 1, MPI_INT, 0 /*master*/, MSG_TAG_WORKER_RESULT, MPI_COMM_WORLD))
    {
        abort_msg("MPI_Send: workerSendResult");
    }
}


std::vector<int> workerReceiveJob() {
    constexpr int tag = MSG_TAG_WORKER_JOB;
    {
        // Non blockin checking if there is message
        int flags = 0;
        while(!flags) {
            if (MPI_SUCCESS != MPI_Iprobe(0 /*master*/, tag, MPI_COMM_WORLD, &flags, MPI_STATUS_IGNORE)){
                abort_msg("workerReceiveJob: MPI_Probe");
            }
            msleep(RECV_FREQ_MS);
        }
    }

    std::vector<int> input_data(MAX_ELEMENTS);
    MPI_Status status;
    if (MPI_SUCCESS != MPI_Recv(input_data.data(), input_data.size(), MPI_INT, 0 /*master*/,
                                tag, MPI_COMM_WORLD, &status)) {
        abort_msg("workerReceiveJob: MPI_Recv");
    }
    int count;
    MPI_Get_count(&status, MPI_INT, &count);
    input_data.resize(count);
    return input_data;
}


class Master
{
    int tile_id = 0;
    int tile_id_max;

    std::vector<int> getTasks(int task_count) {
        std::vector<int> result;
        int i = 0;
        while(i < task_count && tile_id < tile_id_max) {
            result.push_back(tile_id);
            i += 1;
            tile_id += 1;
        }
        return result;
    }

    bool isMoreJobAvailable() {
        return tile_id < tile_id_max;
    }

    void scheduleJobs(int worker_rank, int task_count) {
        auto tasks = getTasks(task_count);
        assert(!tasks.empty());
        masterSendWorkerJob(worker_rank, tasks);
    }

    std::pair<int, int> getNodeResult() {
        constexpr int tag = MSG_TAG_WORKER_RESULT;
        {
            int flags = 0;
            while(!flags) {
                if (MPI_SUCCESS != MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &flags, MPI_STATUS_IGNORE)){
                    abort_msg("workerReceiveJob: MPI_Probe");
                }
                msleep(RECV_FREQ_MS);
            }
        }

        MPI_Status status;
        int tile_id; // tile_id that the worker calculated
        if (MPI_SUCCESS != MPI_Recv(&tile_id, 1, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status)) {
            abort_msg("wait_for_workers: MPI_Send");
        }
        const int input = tile_id;
        const int sender_rank = status.MPI_SOURCE;
        return std::make_pair(sender_rank, input);
    }

public:
    Master(std::vector<WorkerInfo> in_workers, int tasks): tile_id_max{tasks} {
        std::cout << "Master: number of tasks: " << tasks << '\n';
        for (const auto& worker: in_workers) {
            if (workers.find(worker.worker_rank) != workers.cend()) {
                abort_msg("Duplicated worker!");
            }
            workers[worker.worker_rank] = worker;
        }
    }

    void initialSchedule() {
        for (const auto& worker: workers) {
            if (isMoreJobAvailable()) {
                const int tasks = worker.second.cpus * 2.5;
                scheduleJobs(worker.second.worker_rank, tasks);
            } else {
                break;
            }
        }
    }

    void notifyAllWorkerNodesDone() const {
        for (const auto& worker: workers) {
            masterSendWorkerStop(worker.second.worker_rank);
        }
    }

    bool allDone() const {
        return results.size() == tile_id_max;
    }

    void start() {
        initialSchedule();

        while(!allDone()) {
            const auto [sender_rank, input] = getNodeResult();

            printf("Mater: got result job rank: %d x=%d\n", sender_rank, input);
            results.push_back({.worker_rank = sender_rank, .input = input});

            // Send instruction what to do next.
            if (isMoreJobAvailable()) {
                scheduleJobs(sender_rank, 1 /* one task */);
            } else {
                if (!isAllDoneMsgSent) {
                    notifyAllWorkerNodesDone();
                    isAllDoneMsgSent = true;
                }
            }
        }
    }

    bool isAllDoneMsgSent {false};

    std::vector<WorkerResult> results;
    std::map<int, WorkerInfo> workers; // worker_rank -> worker
};



int rayTracerMPI()
{
    // Initialize the MPI environment. The two arguments to MPI Init are not
    // currently used by MPI implementations, but are there in case future
    // implementations might need the arguments.
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    std::vector<WorkerInfo> workers;
    if (world_rank == 0)
    {
        printf("Master node. " macroStringfy(BUILD_TIME) ". Processor %s, rank %d out of %d processors. CPUs: %d   CPUs available: %d\n",
               processor_name, world_rank, world_size, get_nprocs_conf(), get_nprocs());

        // Registration all nodes
        for (int i = 0; i < world_size - 1; i++)
        {
            MPI_Status status;
            int node_info; // CPU count in node
            if (MPI_SUCCESS != MPI_Recv(&node_info, 1, MPI_INT, MPI_ANY_SOURCE, MSG_TAG_REGISTER, MPI_COMM_WORLD, &status))
            {
                abort_msg("MPI_Send");
            }
            WorkerInfo w{status.MPI_SOURCE, node_info}; // rank, CPU count
            workers.push_back(w);
            printf("Registered node: %s\n", w.info().c_str());
        }
    }
    else
    {
        printf("Worker node. " macroStringfy(BUILD_TIME) ". Processor %s, rank %d out of %d processors. CPUs: %d   CPUs available: %d\n",
               processor_name, world_rank, world_size, get_nprocs_conf(), get_nprocs());
        const int nprocs = get_nprocs();
        if (MPI_SUCCESS != MPI_Send(&nprocs, 1, MPI_INT, 0 /*master*/, MSG_TAG_REGISTER, MPI_COMM_WORLD))
        {
            abort_msg("MPI_Send");
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0)
    {
        printf("Master: registration done\n");
        ImageGather image_gather;
        image_gather.gatherTiles();

        Master ms{workers, tiles_columns*tiles_rows};
        ms.start();

        printf("Master results:\n");
        std::map<int, int> stats;
        for(const auto& r: ms.results) {
            stats[r.worker_rank]++;
        }
        for (const auto it: stats) {
            std::cout << "   Node ranks: " <<  it.first << ": " << it.second << '\n';
        }
        printf("Wait for all images parts\n");
        image_gather.join(); // wait until all parts of the image are gathered
        image_gather.image.writeToFileBMP("project/output.bmp");
        printf("Master done\n");
    }
    else
    {
        printf("Worker rank %d: start\n", world_rank);
        ThreadPool pool{get_nprocs()};
        std::atomic_int tasks_count = 0; // for debug purposes (or rather measurement)
        DRenderer renderer;

        //ThreadPool pool_2{1};
        std::mutex mx;

        std::vector<int> inputs;
        bool terminated = false;
        while(!terminated) {
            std::vector<int> input_vector = workerReceiveJob();

            printf("Worker rank %d got tasks: ", world_rank);
            for (int input: input_vector) {
                printf("%d ", input);
            }
            printf("\n");
            for (int input: input_vector) {
                if (input == -1) {
                    terminated = true;
                    break;
                }
                tasks_count++;
                pool.pushTask([input, &renderer, &tasks_count, &mx, &inputs](int thread_id){
                    renderer.renderTile(input);

                    // informa about results
                    workerSendWorkDone(input);
                    //pool.pushTask([input, &renderer](int){
                    //workerSendImageTile(renderer.getImageRef(), input);
                    //});
                    {
                        std::lock_guard lk{mx};
                        inputs.push_back(input);
                    }
                    tasks_count--;
                });
            }

            printf("Worker rank %d. Tasks in buffer: %d\n", world_rank, tasks_count.load());
        }
        // Wait for all threads done
        pool.shutdown();
        printf("Worker rank %d: sending results\n", world_rank);
        for (int input: inputs) {
            workerSendImageTile(renderer.getImageRef(), input);
        }
        printf("Worker rank: %d  done\n", world_rank);
    }

    TileImage::freeMPIType();

    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return 0;
}


int testTileTransfer() {
    // Initialize the MPI environment. The two arguments to MPI Init are not
    // currently used by MPI implementations, but are there in case future
    // implementations might need the arguments.
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    printf(macroStringfy(BUILD_TIME) ". Processor %s, rank %d out of %d processors. CPUs: %d   CPUs available: %d\n",
           processor_name, world_rank, world_size, get_nprocs_conf(), get_nprocs());

    if (world_rank == 0)
    {
        Image image{frame_columns, frame_rows};
        for (int i = 0; i < TileInfo::numberOfTiles(); i++) {
            receiveImageTile(image);
        }
        image.writeToFileBMP("output.bmp");
    }
    else
    {
        Image image{frame_columns, frame_rows};
        test_feelImage(image);
        for (int i = 0; i < TileInfo::numberOfTiles(); i++) {
            workerSendImageTile(image, i);
        }
    }

    //MPI_Barrier(MPI_COMM_WORLD);
    TileImage::freeMPIType();

    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return 0;
}



int main(int, char **) {
    //test_iterateImageTile2();
    //return 0;
    //return testTileTransfer();
    return rayTracerMPI();
}

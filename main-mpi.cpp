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

#define macroStringfy(x) #x

#define MSG_TAG_REGISTER 0
#define MSG_TAG_WORKER_RESULT 1
#define MSG_TAG_WORKER_JOB 2

#define RECV_FREQ_MS 100

constexpr int MAX_ELEMENTS = 300;

void abort_msg(const char msg[])
{
    fprintf(stderr, "ABORT: %s\n", msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}


struct TimeMeasure {
    double full_duration = 0.0; // in seconds
    std::chrono::high_resolution_clock::time_point start_time;

    inline void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    /* Returns the number of seconds elapsed since the last start() call. */
    inline double stop() {
        const auto endv = std::chrono::high_resolution_clock::now();
        if (start_time.time_since_epoch() == std::chrono::high_resolution_clock::time_point::duration::zero()) {
            start_time = std::chrono::high_resolution_clock::time_point(); // zeroing
            return 0.0;
        }
        const auto duration = endv - start_time;
        const double result = duration.count() * 1e-9;
        full_duration += result;
        return result;
    }
};




#include <time.h>
#include <errno.h>

/* msleep(): Sleep for the requested number of milliseconds. */
int msleep(long msec)
{
    struct timespec ts;
    int res;

    if (msec < 0)
    {
        errno = EINVAL;
        return -1;
    }

    ts.tv_sec = msec / 1000;
    ts.tv_nsec = (msec % 1000) * 1000000;

    do {
        res = nanosleep(&ts, &ts);
    } while (res && errno == EINTR);

    return res;
}

std::mutex debug_output;

using Task = std::function<void(int)>;


class ThreadPool {
public:
    ThreadPool(int thread_count) {
        //printf("Worker rank: %d  starts threads: %d\n", world_rank, nprocs);
        for (int i = 0; i < thread_count; i++) {
            threads.push_back(
                std::thread{
                    [this, i](){ workerLoop(i);
                }});
        }
    }

    void pushTask(const Task& t) {
        {
            std::lock_guard<std::mutex> lk{job_mutex};
            job_queue.push_back(t);
        }
        condition.notify_one();
    }

    void shutdown() {
        {
            std::lock_guard<std::mutex> lk{job_mutex};
            terminated = true;
        }
        condition.notify_all(); // wake up all threads
        for(std::thread& th : threads) {
            th.join();
        }
        threads.clear();
    }

    int tasksWaiting() {
         std::lock_guard<std::mutex> lock{job_mutex};
         return job_queue.size();
    }

private:
    void workerLoop(int thread_id) {
        TimeMeasure idle_time;
        while(true) {
            std::unique_lock<std::mutex> lock{job_mutex}; // acquire mutex
            condition.wait(lock, [this]{return !job_queue.empty() || terminated; });
            // We own the lock here

            Task task;
            if (job_queue.empty()) {
                // There is nothing to do
                if (terminated) {
                    break;
                } else {
                    continue;
                }
            } else {
                // There is something to do
                task = job_queue.back();
                job_queue.pop_back();
            }
            lock.unlock();

            idle_time.stop();
            // Execute task
            task(thread_id);
            idle_time.start();
        }

        {
            std::lock_guard<std::mutex> lk{debug_output};
            std::cout << "Thread idle time: " << idle_time.full_duration << " sec\n";
        }
    }

    std::vector<Task> job_queue;
    bool terminated = false;
    std::condition_variable condition;
    std::mutex job_mutex;
    std::vector<std::thread> threads;
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
    int result;
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


void workerSendResult(const int* data) {
    if (MPI_SUCCESS != MPI_Send(data, 2, MPI_INT, 0 /*master*/, MSG_TAG_WORKER_RESULT, MPI_COMM_WORLD))
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

    std::tuple<int, int, int> getNodeResult() {
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
        int result_raw[2];
        if (MPI_SUCCESS != MPI_Recv(&result_raw, 2, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status)) {
            abort_msg("wait_for_workers: MPI_Send");
        }
        const int input = result_raw[0];
        const int result = result_raw[1];
        const int sender_rank = status.MPI_SOURCE;
        return std::make_tuple(sender_rank, input, result);
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
            const auto [sender_rank, input, result] = getNodeResult();

            printf("Mater: got result job rank: %d x=%d y=%d\n", sender_rank, input, result);
            results.push_back({.worker_rank = sender_rank, .input = input, .result = result});

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
        Master ms{workers, DRenderer::numberOfTiles()};
        ms.start();
        printf("Master results:\n");
        std::map<int, int> stats;
        for(const auto& r: ms.results) {
            stats[r.worker_rank]++;
        }
        for (const auto it: stats) {
            std::cout << "Node ranks: " <<  it.first << ": " << it.second << '\n';
        }
        printf("Master done\n");
    }
    else
    {
        printf("Worker rank %d: start\n", world_rank);
        ThreadPool pool{get_nprocs()};
        std::atomic_int tasks_count; // for debug purposes (or rather measurement)
        DRenderer renderer;

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
                pool.pushTask([input, &renderer, &tasks_count](int thread_id){
                      renderer.renderTile(input);

//                    {
//                        std::lock_guard<std::mutex> lk{debug_output};
//                        std::cout << "Thread: " << thread_id << " task: " << input << "  START\n";
//                    }
//                        int r = std::rand() % 4 + 1;
//                        sleep(r);
//                        int y = input * input;
//                    {
//                        std::lock_guard<std::mutex> lk{debug_output};
//                        std::cout << "Thread: " << thread_id << " task: " << input << "->" << y << "  DONE\n";
//                    }

                    // informa about results
                    int tmp[2];
                    tmp[0] = input;
                    tmp[1] = 0;
                    workerSendResult(tmp);
                    tasks_count--;
                });
            }

            printf("Worker rank %d. Tasks in buffer: %d\n", world_rank, tasks_count.load());
        }
        printf("Worker rank %d: stop\n", world_rank);
        // Wait for all threads done
        pool.shutdown();
        renderer.writeImage();
        printf("Worker rank: %d  done\n", world_rank);
    }

    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
    return 0;
}


int main(int, char **) {
    return rayTracerMPI();
}

#pragma once
#include <functional>
#include <mutex>
#include <thread>
#include <condition_variable>
#include "timing.h"

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
            //std::lock_guard<std::mutex> lk{debug_output};
            //std::cout << "Thread idle time: " << idle_time.full_duration << " sec\n";

            // TODO: add gathering time results
        }
    }

    std::vector<Task> job_queue;
    bool terminated = false;
    std::condition_variable condition;
    std::mutex job_mutex;
    std::vector<std::thread> threads;

};

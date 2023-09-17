#include <iostream>
#include <thread>
#include <queue>
#include <functional>
#include <vector>
#include <mutex>
#include <future>
#include <condition_variable>

class ThreadPool {
public:
     ThreadPool(size_t numThreads) {
        for (size_t i = 0; i < numThreads; ++i) {
            threads.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(mutex);

                        // Wait for a task if the queue is empty
                        condition.wait(lock, [this]() { return !tasks.empty(); });

                        // Get the task from the queue
                        task = std::move(tasks.front());
                        tasks.pop();
                    }

                    // Execute the task
                    if (task) {
                         //while (!cancelFlag.load() && !taskCompletedFlag.load()) 
                            task();
                      
                    }
                }
            });
        }
    }


    ~ThreadPool() {
        // Join all threads in the pool

        //cancelFlag.store(true);

        for (auto& thread : threads) {
            thread.join();
        }


    }


    template <class Function, class... Args>
    auto enqueue(Function&& f, Args&&... args)
        -> std::future<decltype(f(std::forward<Args>(args)...))> {
        using return_type = decltype(f(std::forward<Args>(args)...));

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<Function>(f), std::forward<Args>(args)...));

        std::future<return_type> result = task->get_future();

        {
            std::unique_lock<std::mutex> lock(mutex);
            tasks.emplace([task]() { (*task)(); });
        }

        // Notify one thread to execute the task
        condition.notify_one();

        return result;
    }
    //void handleTimeout(std::future<void>& future) {
        // Set the taskCompletedFlag to true when a timeout occurs
     //   taskCompletedFlag.store(true);
    //}

private:
    std::vector<std::thread> threads;
    std::queue<std::function<void()>> tasks;
    std::mutex mutex;
    std::condition_variable condition;
    //std::atomic<bool> cancelFlag{false};
    //std::atomic<bool> taskCompletedFlag{false};
};


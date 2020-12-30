#include "random.h"
#include <random>
#include <chrono>
#include <iostream>

// TODO: This should not be a class anymore but pure module (just set of functions + state)

class RandomGenerator {
public:
    RandomGenerator() {
        //displaySeed();
    }

//    RandomGenerator(int seed): seed{seed} {
//        //displaySeed();
//    }

    void displaySeed() {
        std::cout << "Seed: " << this->seed << '\n';
    }

    float get1() {
        return distribution1(generator);
    }
    float get01() {
        return distribution01(generator);
    }
    int get_int(int i) {
        return dist_int(generator) % i;
    }
private:
    const int seed = std::chrono::system_clock::now().time_since_epoch().count();
    //const long seed = 0; // Seed for testing
    std::default_random_engine generator;
    //static std::uniform_real_distribution<double> distribution05(-0.5f, 0.5f);
    std::uniform_real_distribution<float> distribution1{-1.0, 1.0};
    std::uniform_real_distribution<float> distribution01{0.0, 1.0};
    std::uniform_int_distribution<int> dist_int{0, 1000000}; // TODO: find something better
};


namespace randomns {
thread_local RandomGenerator rnd;

float get1() {
    return rnd.get1();
}
float get01() {
    return rnd.get01();
}
} // namespace: random

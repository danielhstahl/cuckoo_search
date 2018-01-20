#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include <iostream>
#include "cuckoo.h"

TEST_CASE("Test getBestNest 1", "[Cuckoo]"){
    std::vector<std::vector<double> > nest(2, std::vector<double>(2, 1.0));
    std::vector<std::vector<double> > newNest(2, std::vector<double>(2, 4.0));
    std::vector<double> fitness(2, 2.0);
    std::vector<double> best(2, 50000);
    double fmin=cuckoo::getBestNest(&nest, &fitness, &best, newNest, [](const auto& inputs){
        return inputs[0]*inputs[0]+inputs[1]*inputs[1];
    });
    REQUIRE(fitness==std::vector<double>({2.0, 2.0}));
    REQUIRE(best==std::vector<double>({1.0, 1.0}));
    REQUIRE(fmin==2.0);
    REQUIRE(nest==std::vector<std::vector<double> >({
        std::vector<double>({1.0, 1.0}), 
        std::vector<double>({1.0, 1.0})
    }));
}
TEST_CASE("Test getBestNest 2", "[Cuckoo]"){
    std::vector<std::vector<double> > nest(2, std::vector<double>(2, 1.0));
    std::vector<std::vector<double> > newNest(2, std::vector<double>(2, 4.0));
    nest[0][0]=3;
    nest[1][0]=5;
    nest[0][1]=4;
    nest[1][1]=2;
    newNest[0][0]=4;
    newNest[1][0]=2;
    newNest[0][1]=5;
    newNest[1][1]=3;
    std::vector<double> fitness(2, 2.0);
    std::vector<double> best(2, 50000);
    double fmin=cuckoo::getBestNest(&nest, &fitness, &best, newNest, [](const auto& inputs){
        return inputs[0]*inputs[0]+inputs[1]*inputs[1];
    });
    REQUIRE(fitness==std::vector<double>({2.0, 2.0}));
    REQUIRE(best==std::vector<double>({3.0, 4.0}));
    REQUIRE(fmin==2.0);
    REQUIRE(nest==std::vector<std::vector<double> >({
        std::vector<double>({3.0, 4.0}), 
        std::vector<double>({5.0, 2.0})
    }));
}
TEST_CASE("Test cuckoos 1 ", "[Cuckoo]"){
    std::vector<cuckoo::upper_lower<double> > ul;
    std::vector<std::vector<double> > nest(2, std::vector<double>(2, 1.0));
    std::vector<std::vector<double> > newNest(2, std::vector<double>(2, 1.0));
    std::vector<double> best(2, 3.0);
    cuckoo::upper_lower<double> bounds={-1.0, 1.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    cuckoo::getCuckoos(&newNest, nest, best, ul, 1.5, [](){return .5;}, [](){return 1.5;});
    REQUIRE(newNest[0][0]==Approx(.95238));//from matlab ref
    //std::cout<<nest[0][0]<<", "<<nest[0][1]<<", "<<nest[1][0]<<", "<<nest[1][1]<<std::endl;
}
TEST_CASE("Test cuckoos 2", "[Cuckoo]"){
    std::vector<cuckoo::upper_lower<double> > ul;
    std::vector<std::vector<double> > nest(2, std::vector<double>(2, 1.0));
    std::vector<std::vector<double> > newNest(2, std::vector<double>(2, 1.0));
    std::vector<double> best(2, 3.0);
    cuckoo::upper_lower<double> bounds={-1.0, 1.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    cuckoo::getCuckoos(&newNest, nest, best, ul, 1.5, [](){return .5;}, [](){return -1.5;});
    REQUIRE(newNest[0][0]==1.0);//from matlab ref
    //std::cout<<nest[0][0]<<", "<<nest[0][1]<<", "<<nest[1][0]<<", "<<nest[1][1]<<std::endl;
}

TEST_CASE("Test Simple Function", "[Cuckoo]"){
    std::vector<cuckoo::upper_lower<double> > ul;
    cuckoo::upper_lower<double> bounds={-4.0, 4.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
   // ul.push_back(bounds);
    auto results=cuckoo::optimize([](const std::vector<double>& inputs){
        return inputs[0]*inputs[0]+inputs[1]*inputs[1]+inputs[2]*inputs[2]+inputs[3]*inputs[3];
    }, ul, 20, 1000, 42);
    auto params=std::get<cuckoo::optparms>(results);
    std::cout<<params[0]<<", "<<params[1]<<std::endl;
    REQUIRE(std::get<cuckoo::fnval>(results)==Approx(0.0));
}  
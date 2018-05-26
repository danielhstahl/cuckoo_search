#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include <iostream>
#include "firefly.h"
#include "cuckoo.h"

TEST_CASE("Test Simple Function", "[Cuckoo]"){
    std::vector<swarm_utils::upper_lower<double> > ul;
    swarm_utils::upper_lower<double> bounds={-4.0, 4.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    auto results=cuckoo::optimize([](const std::vector<double>& inputs){
        return inputs[0]*inputs[0]+inputs[1]*inputs[1]+inputs[2]*inputs[2]+inputs[3]*inputs[3];
    }, ul, 25, 1000, .00000001, 42);
    auto params=std::get<swarm_utils::optparms>(results);
    std::cout<<params[0]<<", "<<params[1]<<std::endl;
    REQUIRE(std::get<swarm_utils::fnval>(results)==Approx(0.0));
}  
TEST_CASE("Test Rosenbrok Function", "[Cuckoo]"){
    std::vector<swarm_utils::upper_lower<double> > ul;
    swarm_utils::upper_lower<double> bounds={-4.0, 4.0};
    ul.push_back(bounds);
    ul.push_back(bounds);

    auto results=cuckoo::optimize([](const std::vector<double>& inputs){
        return futilities::const_power(1-inputs[0], 2)+100*futilities::const_power(inputs[1]-futilities::const_power(inputs[0], 2), 2);
    }, ul, 20, 10000, .00000001, 42);
    auto params=std::get<swarm_utils::optparms>(results);
    //std::cout<<params[0]<<", "<<params[1]<<std::endl;
    REQUIRE(std::get<swarm_utils::fnval>(results)==Approx(0.0));
    //REQUIRE(params[0]==Approx(1.0));
    //REQUIRE(params[1]==Approx(1.0));
}  
TEST_CASE("Test u^2 Function", "[Cuckoo]"){
    std::vector<swarm_utils::upper_lower<double> > ul;
    swarm_utils::upper_lower<double> bounds={-5.0, 5.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    int numNests=25;
    int maxMC=25000;
    auto results=cuckoo::optimize([](const std::vector<double>& inputs){
        return futilities::sum(inputs, [](const auto& v, const auto& index){
            return futilities::const_power(v-1.0, 2);
        });
    }, ul, numNests, maxMC, .00000001, 42);
    auto params=std::get<swarm_utils::optparms>(results);
    //std::cout<<params[0]<<", "<<params[1]<<std::endl;
    std::cout<<"obj fn: "<<std::get<swarm_utils::fnval>(results)<<std::endl;
    REQUIRE(std::get<swarm_utils::fnval>(results)==Approx(0.0));
    //REQUIRE(params[0]==Approx(1.0));
    //REQUIRE(params[1]==Approx(1.0));
}  
constexpr double rastigrinScale=10;
TEST_CASE("Test Rastigrin Function", "[Cuckoo]"){
    std::vector<swarm_utils::upper_lower<double> > ul;
    swarm_utils::upper_lower<double> bounds={-4.0, 4.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    int n=ul.size();
    auto results=cuckoo::optimize([](const std::vector<double>& inputs){
        return rastigrinScale*inputs.size()+futilities::sum(inputs, [](const auto& val, const auto& index){
            return futilities::const_power(val, 2)-rastigrinScale*cos(2*M_PI*val);
        });
    }, ul, 25, 10000, .00000001, 42);
    auto params=std::get<swarm_utils::optparms>(results);
    for(auto& v:params){
        std::cout<<v<<",";
    }
    REQUIRE(std::get<swarm_utils::fnval>(results)==Approx(0.0));
    //REQUIRE(params[0]==Approx(1.0));
    //REQUIRE(params[1]==Approx(1.0));
}  

/**WOW this works well*/

TEST_CASE("Test Simple Function Firefly", "[FireFly]"){
    std::vector<swarm_utils::upper_lower<double> > ul;
    swarm_utils::upper_lower<double> bounds={-4.0, 4.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    auto results=firefly::optimize([](const std::vector<double>& inputs){
        return inputs[0]*inputs[0]+inputs[1]*inputs[1]+inputs[2]*inputs[2]+inputs[3]*inputs[3];
    }, ul, 1000, 42);
    auto params=std::get<swarm_utils::optparms>(results);
    std::cout<<params[0]<<", "<<params[1]<<std::endl;
    REQUIRE(std::get<swarm_utils::fnval>(results)==Approx(0.0));
}  

/**WOW this works well*/
TEST_CASE("Test Rosenbrok Function FireFly", "[FireFly]"){
    std::vector<swarm_utils::upper_lower<double> > ul;
    swarm_utils::upper_lower<double> bounds={-4.0, 4.0};
    ul.push_back(bounds);
    ul.push_back(bounds);

    auto results=firefly::optimize([](const std::vector<double>& inputs){
        return futilities::const_power(1-inputs[0], 2)+100*futilities::const_power(inputs[1]-futilities::const_power(inputs[0], 2), 2);
    }, ul,  1000, 42);
    auto params=std::get<swarm_utils::optparms>(results);
    //std::cout<<params[0]<<", "<<params[1]<<std::endl;
    REQUIRE(std::get<swarm_utils::fnval>(results)==Approx(0.0));
    //REQUIRE(params[0]==Approx(1.0));
    //REQUIRE(params[1]==Approx(1.0));
}  
/**WOW this works horribly*/
TEST_CASE("Test u^2 Function for Firefly", "[Firefly]"){
    std::vector<swarm_utils::upper_lower<double> > ul;
    swarm_utils::upper_lower<double> bounds={-5.0, 5.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    int numNests=25;
    int maxMC=1000;
    auto results=firefly::optimize([](const std::vector<double>& inputs){
        return futilities::sum(inputs, [](const auto& v, const auto& index){
            return futilities::const_power(v-1.0, 2);
        });
    }, ul, maxMC, 42);
    auto params=std::get<swarm_utils::optparms>(results);
    std::cout<<"Firefly u^2:"<<std::endl;
    for(auto& v:params){
        std::cout<<v<<",";
    }
    std::cout<<"obj fn: "<<std::get<swarm_utils::fnval>(results)<<std::endl;
    REQUIRE(std::get<swarm_utils::fnval>(results)==Approx(0.0));
    //REQUIRE(params[0]==Approx(1.0));
    //REQUIRE(params[1]==Approx(1.0));
}  

TEST_CASE("Test Rastigrin Function FireFly", "[FireFly]"){
    std::vector<swarm_utils::upper_lower<double> > ul;
    swarm_utils::upper_lower<double> bounds={-4.0, 4.0};
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    ul.push_back(bounds);
    int n=ul.size();
    auto results=firefly::optimize([](const std::vector<double>& inputs){
        return rastigrinScale*inputs.size()+futilities::sum(inputs, [](const auto& val, const auto& index){
            return futilities::const_power(val, 2)-rastigrinScale*cos(2*M_PI*val);
        });
    }, ul, 1000, 42);

    auto params=std::get<swarm_utils::optparms>(results);
    std::cout<<"Firefly Rastigrin:"<<std::endl;
    for(auto& v:params){
        std::cout<<v<<",";
    }
    REQUIRE(std::get<swarm_utils::fnval>(results)==Approx(0.0));
    //REQUIRE(params[0]==Approx(1.0));
    //REQUIRE(params[1]==Approx(1.0));
}  
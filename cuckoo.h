#ifndef __CUCKOO__H__
#define __CUCKOO__H__
#include "FunctionalUtilities.h"
#include <cstdlib> 
#include "SimulateNorm.h"
#include <tuple>

/**Based off the following paper: http://www.airccse.org/journal/ijaia/papers/0711ijaia04.pdf*/
namespace cuckoo{
    constexpr int optparms=0;
    constexpr int fnval=1;
    template<typename T>
    struct upper_lower{
        T upper;
        T lower;
        upper_lower(T&& lower_, T&& upper_){
            upper=upper_;
            lower=lower_;
        }
    };

    template<typename T, typename U>
    auto getRandomParameter(const T& lower, const T& upper, const U& rand){
        return lower+rand*(upper-lower);
    }

    template<typename T, typename U>
    auto getTruncatedParameter(const T& lower, const T& upper, const U& result){
        return result>upper?upper:(result<lower?lower:result);
    }

    template<typename T,  typename U>
    auto getLevy(const T& alpha, const U& rand){
        return pow(rand, -1.0/alpha);
    }

    template<typename T, typename U>
    auto getLevyFlight(const T& currVal, const T& stepSize, const T& lambda, U&& rand, U&& normRand){
        return currVal+stepSize*getLevy(lambda, rand)*normRand;
    }

    auto getUniform(){
        return (double)rand()/RAND_MAX;
    }
    auto getDiscreteUniform(int n){
        return (int)getUniform()*n;
    }
    template<typename Array>
    auto getRandomParameters(const Array& ul){
        return futilities::for_each(0, (int)ul.size(), [&](const auto& index){
            return getRandomParameter(ul[index].lower, ul[index].upper, getUniform());
        });
    }
    template<typename Array, typename Parms, typename U, typename Unif, typename Norm>
    auto getRandomLevyParameters(
        const Array& ul, const Parms& parameters, 
        const U& lambda, U&& alpha,  
        const Unif& unif, const Norm& norm
    ){
        return futilities::for_each(0, (int)ul.size(), [&](const auto& index){
            return getTruncatedParameter(
                ul[index].lower, ul[index].upper, 
                getLevyFlight(
                    parameters[index], alpha, lambda, unif(), norm()
                )
            );
        });
    }

    template<typename Nest>
    void sortNest(Nest& nestRef){
        std::sort(nestRef.begin(), nestRef.end(), [](const auto& val1, const auto& val2){
            return val1.second<val2.second;//smallest to largest
        });
    }
    template<typename Array, typename ObjFn, typename RandParams>
    auto getNewParameterAndFn(const Array& ul, const ObjFn& objFn, const RandParams& rndParms){
        auto parameters=rndParms(ul);
        return std::pair<std::vector<double>, double>(parameters, objFn(parameters));
    }

    template<
        typename Nest,  typename Array, typename ObjFun,
        typename U, typename Norm, typename Unif
    >
    void getCuckoos(
        Nest* nest, 
        const ObjFun& objFun,
        const Array& ul, 
        const U& lambda, 
        U&& alpha, 
        const Unif& unif,
        const Norm& norm
    ){
        Nest& nestRef= *nest;
        int n=nestRef.size(); //number of nests
        int m=nestRef[0].first.size(); //number of parameters
        
        auto i=getDiscreteUniform(n);//get random parameter set
        auto getLevyParameters=[&](const auto& ult){
            return getRandomLevyParameters(ult, nestRef[i].first, lambda, std::move(alpha), unif, norm);
        };
        auto levyParameters=getNewParameterAndFn(ul, objFun, getLevyParameters);
        auto j=getDiscreteUniform(n);//get random nest
        if(nestRef[j].second>levyParameters.second){
            nestRef[j]=levyParameters;
        }
    }


    
    template<typename Array, typename ObjFn, typename GetUniformParameters>
    auto getNewNest(const Array& ul, const ObjFn& objFn, int n, const GetUniformParameters& getUniformParameters){
        return futilities::for_each(0, n, [&](const auto& index){
            return getNewParameterAndFn(ul, objFn, getUniformParameters);
        });
    }

    template<typename P, typename Index>
    P getPA(const P& pMin, const P& pMax, const Index& index, const Index& n){
        return pMax-(pMax-pMin)*index/n;
    }

    template<typename Nest, typename ObjFn, typename P, typename Array, typename GetUniformParameters>
    void emptyNests(
        Nest* newNest, const ObjFn& objFn, const Array& ul, 
        const P& p, const GetUniformParameters& getUniformParameters
    ){
        Nest& nestRef= *newNest;
        int n=nestRef.size();
        int numToKeep=(int)(p*nestRef.size());
        int startNum=n-numToKeep;
        
        for(int i=startNum; i<n; ++i){
            nestRef[i]=getNewParameterAndFn(ul, objFn, getUniformParameters);
        }
    }

    double getAlpha(int index, int totalMC, double alphaMin, double alphaMax){
        double c=log(alphaMin/alphaMax)/totalMC;
        return alphaMax*exp(c*index);
    }


    template< typename Array, typename ObjFn>
    auto optimize(const ObjFn& objFn, const Array& ul, int n, int totalMC, double tol, int seed){
        int numParams=ul.size();
        srand(seed);
        auto getUniformParameters=[](const auto& ult){
            return getRandomParameters(ult);
        };
        auto nest=getNewNest(ul, objFn, n, getUniformParameters);
        //auto newNest=getNewNest(ul, objFn, n);
        double lambda=1.5;
        double alphaMin=.01;
        double alphaMax=.5;
        double pMin=.05;
        double pMax=.5;
        SimulateNorm norm(seed);
        
        double fMin=2;
        int i=0;
        auto unifL=[](){return getUniform();};
        auto normL=[&](){return norm.getNorm();};
        while(i<totalMC&&fMin>tol){
            /**Completely overwrites newNest*/
            getCuckoos(&nest, objFn, ul, 
                lambda, 
                getAlpha(i, totalMC, alphaMin, alphaMax),
                unifL, 
                normL
            );
            sortNest(nest);
            emptyNests(&nest, objFn, ul, getPA(pMin, pMax, i, totalMC), getUniformParameters);
            sortNest(nest);
            fMin=nest[0].second;

            #ifdef VERBOSE_FLAG
                std::cout<<"Index: "<<i<<", Param Vals: ";
                for(auto& v:nest[0].first){
                    std::cout<<v<<", ";
                }
                std::cout<<", Obj Val: "<<fMin<<std::endl;
            #endif
            ++i;
        }
        return nest[0];
    }
}



#endif
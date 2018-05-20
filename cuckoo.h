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
    template<typename Array>
    auto getRandomParameters(const Array& ul){
        return futilities::for_each(0, (int)ul.size(), [&](const auto& index){
            return getRandomParameter(ul[index].lower, ul[index].upper, getUniform());
        });
    }

    template<typename Nest>
    void sortNest(Nest& nestRef){
        std::sort(nestRef.begin(), nestRef.end(), [](const auto& val1, const auto& val2){
            return val1.second<val2.second;//smallest to largest
        });
    }
    template<typename Nest>
    void getBestNest(Nest* nest, const Nest& newNest){
        Nest& nestRef= *nest;
        for(int i=0; i<nestRef.size(); ++i){
            if(newNest[i].second<=nestRef[i].second){
                nestRef[i].second=newNest[i].second;//replace previous function result with current
                nestRef[i].first=newNest[i].first; //replace previous parameters with current
            }
        }
        sortNest(nestRef);
 
    }

    template<
        typename Nest,  typename Array, typename ObjFun,
        typename U, typename Norm, typename Unif
    >
    void getCuckoos(
        Nest* newNest, const Nest& nest, 
        const ObjFun& objFun,
        const Array& ul, 
        const U& lambda, 
        U&& alpha, 
        const Unif& unif,
        const Norm& norm
    ){
        int n=nest.size();
        int m=nest[0].first.size();
        Nest& nestRef= *newNest;
        for(int i=0; i<n;++i){
            for(int j=0; j<m; ++j){
                nestRef[i].first[j]=getTruncatedParameter(
                    ul[j].lower, ul[j].upper, 
                    getLevyFlight(
                        nest[i].first[j], alpha, lambda, unif(), norm()
                    )
                );
            }
            nestRef[i].second=objFun(nestRef[i].first);
        }
    }


    template<typename Array, typename ObjFn>
    auto getNewParameterAndFn(const Array& ul, const ObjFn& objFn){
        auto parameters=getRandomParameters(ul);
        return std::pair<std::vector<double>, double>(parameters, objFn(parameters));
    }
    template<typename Array, typename ObjFn>
    auto getNewNest(const Array& ul, const ObjFn& objFn, int n){
        return futilities::for_each(0, n, [&](const auto& index){
            return getNewParameterAndFn(ul, objFn);
        });
    }

    template<typename P, typename Index>
    P getPA(const P& pMin, const P& pMax, const Index& index, const Index& n){
        return pMax-(pMax-pMin)*index/n;
    }

    template<typename Nest, typename ObjFn, typename P, typename Array>
    void emptyNests(Nest* newNest, const ObjFn& objFn, const Array& ul, const P& p){
        Nest& nestRef= *newNest;
        int n=nestRef.size();
        int numToKeep=(int)(p*nestRef.size());
        int startNum=n-numToKeep;
        for(int i=startNum; i<n; ++i){
            nestRef[i]=getNewParameterAndFn(ul, objFn);
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
        auto nest=getNewNest(ul, objFn, n);
        auto newNest=getNewNest(ul, objFn, n);
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
            //newNest now has the previous values from nest with levy flights added
            getCuckoos(&newNest, nest, objFn, ul, 
                lambda, 
                getAlpha(i, totalMC, alphaMin, alphaMax),
                unifL, 
                normL
            );
            //compare previous nests with cuckoo nests and sort results
            //nest now has the best of nest and newNest
            getBestNest(
                &nest, 
                newNest
            );
            //remove bottom "p" nests and resimulate.
            emptyNests(&nest, objFn, ul, getPA(pMin, pMax, i, totalMC));
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
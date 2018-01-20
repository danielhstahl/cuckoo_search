#ifndef __CUCKOO__H__
#define __CUCKOO__H__
#include "FunctionalUtilities.h"
#include <cstdlib> 
#include "SimulateNorm.h"

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
    auto getLevyFlight(const T& currVal, const T& currBest, const T& stepSize, const T& alpha, U&& rand, U&& normRand){
        return currVal+stepSize*(currVal-currBest)*getLevy(alpha, rand)*normRand;
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
    template<typename Nest, typename ObjFn, typename FnResultArray, typename BestParameter>
    double getBestNest(Nest* nest, FnResultArray* fnResult, BestParameter* bestParam, const Nest& newNest, const ObjFn& objFn){
        Nest& nestRef= *nest;
        FnResultArray& fnResultRef= *fnResult;
        BestParameter& bestParamRef= *bestParam;
        for(int i=0; i<nestRef.size(); ++i){
            auto newFnResult=objFn(newNest[i]);
            if(newFnResult<=fnResultRef[i]){
                fnResultRef[i]=newFnResult;
                nestRef[i]=newNest[i];
            }
        }
        double min=1000000;
        int indexOfMin=0;
        for(int i=0; i<nestRef.size(); ++i){
            if(min>fnResultRef[i]){
                min=fnResultRef[i];
                indexOfMin=i;
            }
        }
        bestParamRef=nestRef[indexOfMin];
        return min;
    }
    template<typename Nest, typename ObjFn, typename FnResultArray, typename BestParameter>
    double getInitBestNest(const Nest& nest, FnResultArray* fnResult, BestParameter* bestParam, const ObjFn& objFn){
        FnResultArray& fnResultRef= *fnResult;
        BestParameter& bestParamRef= *bestParam;
        for(int i=0; i<nest.size(); ++i){
            auto newFnResult=objFn(nest[i]);
            if(newFnResult<=fnResultRef[i]){
                fnResultRef[i]=newFnResult;
            }
        }
        double min=1000000;
        int indexOfMin=0;
        for(int i=0; i<nest.size(); ++i){
            if(min>fnResultRef[i]){
                min=fnResultRef[i];
                indexOfMin=i;
            }
        }
        bestParamRef=nest[indexOfMin];
        return min;
    }

    template<
        typename Nest, typename BestParameter, typename Array, 
        typename U, typename Norm, typename Unif
    >
    void getCuckoos(
        Nest* newNest, const Nest& nest, const BestParameter& best, 
        const Array& ul, 
        const U& alpha, 
        Unif&& unif,
        Norm&& norm
    ){
        int n=nest.size();
        Nest& nestRef= *newNest;
        for(int i=0; i<n;++i){
            for(int j=0; j<best.size(); ++j){
                nestRef[i][j]=getTruncatedParameter(
                    ul[j].lower, ul[j].upper, 
                    getLevyFlight(
                        nest[i][j], best[j], .01, alpha, unif(), norm()
                    )
                );
            }
        }
    }


    template<typename Array>
    auto getNewNest(const Array& ul, int n){
        return futilities::for_each(0, n, [&](const auto& index){
            return getRandomParameters(ul);
        });
    }

    template<typename P, typename Index>
    P getPA(const P& pMin, const P& pMax, const Index& index, const Index& n){
        return pMax-(pMax-pMin)*index/n;
    }

    template<typename Nest, typename P, typename Array>
    void emptyNests(Nest* newNest, const Nest& nest, const Array& ul, const P& p){
        int n=nest.size();
        Nest& nestRef= *newNest;
        for(int i=0; i<n; ++i){
            if(getUniform()>p){
                nestRef[i]=getRandomParameters(ul);
            }
        }
    }


    template< typename Array, typename ObjFn>
    auto optimize(const ObjFn& objFn, const Array& ul, int n, int totalMC, int seed){
        int numParams=ul.size();
        srand(seed);
        auto nest=getNewNest(ul, n);
        auto newNest=getNewNest(ul, n);
        double alpha=1.5;
        double pMin=.05;
        double pMax=.95;
        SimulateNorm norm;
        std::vector<double> fnResult(n, 10000000);
        std::vector<double> bestParams(numParams, 10000000);
        std::vector<double> tmpParams(numParams, 10000000);

        double fMin=getInitBestNest(nest, &fnResult, &bestParams, objFn);
        std::cout<<fnResult[0]<<", "<<fnResult[n-1]<<std::endl;
        std::cout<<fMin<<std::endl;
        std::cout<<bestParams[0]<<", "<<bestParams[1]<<std::endl;
        double fNew;
        for(int i=0; i<totalMC; ++i){
            /**Completely overwrites newNest*/
            getCuckoos(&newNest, nest, bestParams, ul, alpha, 
                [](){return getUniform();}, 
                [&](){return norm.getNorm();}
            );
            fNew=getBestNest(
                &nest, 
                &fnResult, 
                &tmpParams,//not used until 189...get overwritten below but thats fine
                newNest,
                objFn
            );
            emptyNests(&newNest, nest, ul, getPA(pMin, pMax, i, totalMC));
            fNew=getBestNest(
                &nest, 
                &fnResult, 
                &tmpParams,
                newNest,
                objFn
            );
            if(fNew<fMin){
                fMin=fNew;
                bestParams=tmpParams;
            }
        }
        return std::make_tuple(bestParams, fMin);
    }
}



#endif
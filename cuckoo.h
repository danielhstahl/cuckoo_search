#ifndef __CUCKOO__H__
#define __CUCKOO__H__
#include "FunctionalUtilities.h"
#include <cstdlib> 
#include "SimulateNorm.h"
#include <tuple>
#include "utils.h"

/**Based off the following paper: http://www.airccse.org/journal/ijaia/papers/0711ijaia04.pdf*/
namespace cuckoo{
    
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
    template<typename Parm>
    auto getStepSize(const Parm& curr, const Parm& best, const Parm& lower, const Parm& upper){
        return .01*(upper-lower)*(curr-best);//.01 comes from matlab code
    }

    template<
        typename Nest,  typename Array, typename ObjFun,
        typename U, typename Norm, typename Unif,
        typename BestParameter
    >
    void getCuckoos(
        Nest* newNest, const Nest& nest, 
        const BestParameter& bP,
        const ObjFun& objFun,
        const Array& ul, 
        const U& lambda, 
        const Unif& unif,
        const Norm& norm
    ){
        int n=nest.size(); //num nests
        int m=nest[0].first.size(); //num parameters
        Nest& nestRef= *newNest;
        for(int i=0; i<n;++i){
            for(int j=0; j<m; ++j){
                nestRef[i].first[j]=swarm_utils::getTruncatedParameter(
                    ul[j].lower, ul[j].upper, 
                    swarm_utils::getLevyFlight(
                        nest[i].first[j], 
                        getStepSize(nest[i].first[j], bP[j], ul[j].lower, ul[j].upper), 
                        lambda, unif(), norm()
                    )
                );
            }
            nestRef[i].second=objFun(nestRef[i].first);
        }
    }


    
    template<typename Array, typename ObjFn, typename Rand>
    auto getNewNest(const Array& ul, const ObjFn& objFn, const Rand& rnd, int n){
        return futilities::for_each(0, n, [&](const auto& index){
            return swarm_utils::getNewParameterAndFn(ul, objFn, rnd);
        });
    }

    template<typename P, typename Index>
    P getPA(const P& pMin, const P& pMax, const Index& index, const Index& n){
        return pMax-(pMax-pMin)*index/n;
    }

    template<typename Nest, typename ObjFn, typename P, typename Array, typename Rand>
    void emptyNests(Nest* newNest, const ObjFn& objFn, const Rand& rnd, const Array& ul, const P& p){
        Nest& nestRef= *newNest;
        int n=nestRef.size();
        int numToKeep=(int)(p*nestRef.size());
        int startNum=n-numToKeep;
        for(int i=startNum; i<n; ++i){
            nestRef[i]=swarm_utils::getNewParameterAndFn(ul, objFn, rnd);
        }
    }

    template< typename Array, typename ObjFn>
    auto optimize(const ObjFn& objFn, const Array& ul, int n, int totalMC, double tol, int seed){
        int numParams=ul.size();
        srand(seed);
        SimulateNorm norm(seed);
        auto unifL=[](){return swarm_utils::getUniform();};
        auto normL=[&](){return norm.getNorm();};
        auto nest=getNewNest(ul, objFn,normL, n);
        double lambda=1.5;
        double pMin=.05;
        double pMax=.5;
        
        double fMin=2;
        int i=0;

        sortNest(nest);
        auto newNest=getNewNest(ul, objFn,normL, n);
       
       
        while(i<totalMC&&fMin>tol){
            /**Completely overwrites newNest*/
            //newNest now has the previous values from nest with levy flights added
            getCuckoos(
                &newNest, 
                nest, nest[0].first, //the current best nest
                objFn, ul, 
                lambda, 
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
            emptyNests(&nest, objFn, normL, ul, getPA(pMin, pMax, i, totalMC));
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
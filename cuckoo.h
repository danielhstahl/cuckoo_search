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
    /*template<typename Nest>
    void getBestNest(Nest* nest, const Nest& newNest){
        Nest& nestRef= *nest;
        for(int i=0; i<nestRef.size(); ++i){
            if(newNest[i].second<=nestRef[i].second){
                nestRef[i].second=newNest[i].second;//replace previous function result with current
                nestRef[i].first=newNest[i].first; //replace previous parameters with current
            }
        }
        sortNest(nestRef);
 
    }*/

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
            return getRandomLevyParameters(ult, nestRef[i].first, lambda, alpha, unif, norm);
        };
        auto levyParameters=getNewParameterAndFn(ul, objFun, getLevyParameters);
        auto j=getDiscreteUniform(n);//get random nest
        if(nestRef[j].second>levyParameters.second){
            nestRef[j]=levyParameters;
        }
        /*for(int i=0; i<n;++i){
            for(int j=0; j<m; ++j){
                nestRef[i].first[j]=getTruncatedParameter(
                    ul[j].lower, ul[j].upper, 
                    getLevyFlight(
                        nest[i].first[j], alpha, lambda, unif(), norm()
                    )
                );
            }
            nestRef[i].second=objFun(nestRef[i].first);
        }*/
    }


    
    template<typename Array, typename ObjFn>
    auto getNewNest(const Array& ul, const ObjFn& objFn, int n){
        auto getUniformParameters=[](const auto& ult){  //todo!  once source for this function (put in "main")
            return getRandomParameters(ult);
        };
        return futilities::for_each(0, n, [&](const auto& index){
            return getNewParameterAndFn(ul, objFn, getUniformParameters);
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
        auto getUniformParameters=[](const auto& ult){ //todo!  once source for this function (put in "main")
            return getRandomParameters(ult);
        };
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
        auto nest=getNewNest(ul, objFn, n);
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
            //newNest now has the previous values from nest with levy flights added
            getCuckoos(&nest, objFn, ul, 
                lambda, 
                getAlpha(i, totalMC, alphaMin, alphaMax),
                unifL, 
                normL
            );
            //compare previous nests with cuckoo nests and sort results
            //nest now has the best of nest and newNest
            /*getBestNest(
                &nest, 
                newNest
            );*/
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
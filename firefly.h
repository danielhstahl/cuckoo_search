#ifndef __FIREFLY_H__
#define __FIREFLY_H__
#include "FunctionalUtilities.h"
#include "utils.h"
#include "SimulateNorm.h"
namespace firefly{
    constexpr double beta=1;
    constexpr int n=25;
    template<typename FireFlies>
    void sortNest(FireFlies& fireflyRef){
        std::sort(fireflyRef.begin(), fireflyRef.end(), [](const auto& val1, const auto& val2){
            return val1.second<val2.second;//smallest to largest
        });
    }
    template<typename Params>
    auto getDistanceSq(const Params& params1, const Params& params2){
        return futilities::sum(params1, [&](const auto& v, const auto& i){
            return futilities::const_power(v-params2[i], 2);
        });
    }
    template<typename Xi, typename Xj>
    auto getNextDetStep(const Xi& xi, const Xj& xj, double step){
        return xi+step*(xj-xi);
    }
    
    template<typename FireFlies, typename Norm, typename ObjFn, typename Array>
    void getUpdate(FireFlies* fireflies, const ObjFn& objFun, const Array& ul, double beta, double gamma, double vol, const Norm& norm){
        FireFlies& firefliesRef= *fireflies;
        const int numFlies=firefliesRef.size(); //num flies
        const int numParams=firefliesRef[0].first.size(); //num parameters
        //fireflies are sorted
        for(int i=0; i<numFlies; ++i){
            for(int j=0; j<numFlies; ++j){
                //minimizing, hence the opposite sign
                if(firefliesRef[j].second<firefliesRef[i].second){
                    const double r=getDistanceSq(firefliesRef[i].first, firefliesRef[j].first);
                    for(int k=0; k<numParams; ++k){
                        firefliesRef[i].first[k]=swarm_utils::getTruncatedParameter(
                            ul[k].lower, ul[k].upper,
                            getNextDetStep(
                                firefliesRef[i].first[k],
                                firefliesRef[j].first[k],beta*exp(-gamma*r)
                            )+vol*norm()*(ul[k].upper-ul[k].lower) //should this be scaled by size of input range?
                        );
                    }
                    firefliesRef[i].second=objFun(firefliesRef[i].first);
                }
            }
        }
    }


    template<typename Array, typename ObjFn, typename Rand>
    auto getInitialFirefly(const Array& ul, const ObjFn& objFn, const Rand& rnd, int n){
        return futilities::for_each(0, n, [&](const auto& index){
            return swarm_utils::getNewParameterAndFn(ul, objFn, rnd);
        });
    }


    template< typename Array, typename ObjFn>
    auto optimize(
        const ObjFn& objFn, 
        const Array& ul, 
        int totalMC,  
        int seed
    ){
        srand(seed);
        const int numParams=ul.size();
        const double L=futilities::sum(ul, [](const auto& v, const auto& index){
            return v.upper-v.lower;
        }); //average scale
        const double alpha0=.25;//L*.1;
        const double gamma=1.0/sqrt(L);   
        const double delta=.97;
        double deltaT=delta;
        SimulateNorm norm(seed);
       
        auto unifL=[](){return 2*swarm_utils::getUniform()-1;}; //to keep uniform

        auto fireflies=getInitialFirefly(ul, objFn, unifL, n);
        auto normL=[&](){return norm.getNorm();};
        sortNest(fireflies);
        for(int i=0; i<totalMC; ++i){
            
            getUpdate(&fireflies, objFn, ul, beta, gamma, alpha0*deltaT, normL);
            sortNest(fireflies);
            deltaT*=delta;
            #ifdef VERBOSE_FLAG
                std::cout<<"Index: "<<i<<", Param Vals: ";
                for(auto& v:fireflies[0].first){
                    std::cout<<v<<", ";
                }
                std::cout<<", Obj Val: "<<fireflies[0].second<<std::endl;
            #endif
        }
        return fireflies[0];



    }

}


#endif
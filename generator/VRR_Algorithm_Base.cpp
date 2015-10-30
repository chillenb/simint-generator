#include <iostream>
#include <algorithm>

#include "generator/Helpers.hpp"
#include "generator/VRR_Algorithm_Base.hpp"

using namespace std;

VRR_Algorithm_Base::VRR_Algorithm_Base(const OptionsMap & options)
    : options_(options)
{
}

int VRR_Algorithm_Base::GetOption(int opt) const
{
    return options_.at(opt);
}

QAMList VRR_Algorithm_Base::GetAMOrder(void) const
{
    return amorder_; 
}

QAMSet VRR_Algorithm_Base::GetAllAM(void) const
{
    return allam_;
}

int VRR_Algorithm_Base::GetMReq(QAM am) const
{
    return vrrmreq_max_.at(am);
}

VRR_StepList VRR_Algorithm_Base::GetSteps(QAM am) const
{
    return vrrmap_.at(am);
}


QAMSet VRR_Algorithm_Base::GetAMReq(QAM am) const
{
    return qamreq_.at(am);
}


IntSet VRR_Algorithm_Base::GetIntReq_2p(QAM am) const
{
    if(qamint_2p_.count(am) == 0)
        return IntSet();
    else
        return qamint_2p_.at(am);
}


IntSet VRR_Algorithm_Base::GetIntReq_2q(QAM am) const
{
    if(qamint_2q_.count(am) == 0)
        return IntSet();
    else
        return qamint_2q_.at(am);
}


IntSet VRR_Algorithm_Base::GetIntReq_2pq(QAM am) const
{
    if(qamint_2pq_.count(am) == 0)
        return IntSet();
    else
        return qamint_2pq_.at(am);
}

IntSet VRR_Algorithm_Base::GetAllInt_2p(void) const
{
    IntSet iset;
    for(const auto & it : qamint_2p_)
        iset.insert(it.second.begin(), it.second.end());
    return iset;
}

IntSet VRR_Algorithm_Base::GetAllInt_2q(void) const
{
    IntSet iset;
    for(const auto & it : qamint_2q_)
        iset.insert(it.second.begin(), it.second.end());
    return iset;
}

IntSet VRR_Algorithm_Base::GetAllInt_2pq(void) const
{
    IntSet iset;
    for(const auto & it : qamint_2pq_)
        iset.insert(it.second.begin(), it.second.end());
    return iset;
}


StringSet VRR_Algorithm_Base::GetVarReq(QAM am) const
{
    return varreq_.at(am);
}


StringSet VRR_Algorithm_Base::GetAllVarReq(void) const
{
    StringSet sset;
    for(const auto & it : varreq_)
        sset.insert(it.second.begin(), it.second.end());
    return sset;
}



int VRR_Algorithm_Base::GetMaxInt(void) const
{
    return maxint_;
}



void VRR_Algorithm_Base::PruneQuartets_(QuartetSet & q) const
{
    QuartetSet qnew;
 
    for(const auto & it : q)
    {
        if(it && it.am() > 0)
            qnew.insert(it);
    }

    q = qnew;
}

bool VRR_Algorithm_Base::HasBraVRR(void) const
{
    for(const auto & it : vrrmap_)
    {
        for(const auto & it2 : it.second)
        {
            if(it2.type == RRStepType::I || it2.type == RRStepType::J)
                return true;
        }
    }

    return false;
}

bool VRR_Algorithm_Base::HasKetVRR(void) const
{
    for(const auto & it : vrrmap_)
    {
        for(const auto & it2 : it.second)
        {
            if(it2.type == RRStepType::K || it2.type == RRStepType::L)
                return true;
        }
    }

    return false;
}

void VRR_Algorithm_Base::AMOrder_AddWithDependencies_(QAMList & order, QAM am) const
{
    // skip if it was already done somewhere
    if(std::find(order.begin(), order.end(), am) != order.end())
        return;

    // skip if it's not done by ET
    if(allam_.count(am) == 0 || am == QAM{0,0,0,0})
        return;

    // get requirements
    QAMSet req = GetAMReq(am);

    for(const auto & it : req)
        AMOrder_AddWithDependencies_(order, it);

    order.push_back(am);
}


void VRR_Algorithm_Base::Create(QAM q)
{
    Create(GenerateInitialQuartetTargets(q));
}



void VRR_Algorithm_Base::Create(const QuartetSet & q)
{
    // holds the requirements for each am
    // so store what we initially want
    PrintQuartetSet(q, "Initial VRR");

    VRR_StepMap newmap;

    // solved quartets
    QuartetSet solvedquartets;

    // unsolved quartets
    QuartetSet targets = q;
    PruneQuartets_(targets);


    // add max m for initial targets
    // (should be 0)
    for(const auto & it : q)
    {
        QAM qam = it.amlist();
        if( (vrrmreq_max_.count(qam) == 0) || (it.m > vrrmreq_max_.at(qam)))
            vrrmreq_max_[qam] = it.m;
    }


    while(targets.size())
    {
        QuartetSet newtargets;

        for(const auto & it : targets)
        {
            VRRStep vs = VRRStep_(it);
            QAM qam = it.amlist();

            // set the min/max m value
            if((vrrmreq_max_.count(qam) == 0) || (it.m > vrrmreq_max_.at(qam)))
                vrrmreq_max_[qam] = it.m;
            if((vrrmreq_min_.count(qam) == 0) || (it.m < vrrmreq_min_.at(qam)))
                vrrmreq_min_[qam] = it.m;

            // fill in the ijkl members
            int istep = XYZStepToIdx(vs.xyz);
            vs.ijkl = {0, 0, 0, 0};
            if(vs.type == RRStepType::I)
            {
                if(vs.src[2] || vs.src[3])
                    vs.ijkl[0] = it.bra.left.ijk[istep]-1;
                if(vs.src[4] || vs.src[5])
                    vs.ijkl[1] = it.bra.right.ijk[istep];
                if(vs.src[6])
                    vs.ijkl[2] = it.ket.left.ijk[istep];
                if(vs.src[7])
                    vs.ijkl[3] = it.ket.right.ijk[istep];
            }
            else if(vs.type == RRStepType::J)
            {
                if(vs.src[2] || vs.src[3])
                    vs.ijkl[0] = it.bra.left.ijk[istep];
                if(vs.src[4] || vs.src[5])
                    vs.ijkl[1] = it.bra.right.ijk[istep]-1;
                if(vs.src[6])
                    vs.ijkl[2] = it.ket.left.ijk[istep];
                if(vs.src[7])
                    vs.ijkl[3] = it.ket.right.ijk[istep];
            }
            else if(vs.type == RRStepType::K)
            {
                if(vs.src[2] || vs.src[3])
                    vs.ijkl[2] = it.ket.left.ijk[istep]-1;
                if(vs.src[4] || vs.src[5])
                    vs.ijkl[3] = it.ket.right.ijk[istep];
                if(vs.src[6])
                    vs.ijkl[0] = it.bra.left.ijk[istep];
                if(vs.src[7])
                    vs.ijkl[1] = it.bra.right.ijk[istep];
            }
            else
            {
                if(vs.src[2] || vs.src[3])
                    vs.ijkl[2] = it.ket.left.ijk[istep];
                if(vs.src[4] || vs.src[5])
                    vs.ijkl[3] = it.ket.right.ijk[istep]-1;
                if(vs.src[6])
                    vs.ijkl[0] = it.bra.left.ijk[istep];
                if(vs.src[7])
                    vs.ijkl[1] = it.bra.right.ijk[istep];
            }
            
            // add new targets and m values 
            for(const auto & it2 : vs.src)
            {
                if(it2)
                {
                    if(solvedquartets.count(it2) == 0)
                        newtargets.insert(it2);

                    QAM qam2 = it2.amlist();
                    if( (vrrmreq_max_.count(qam2) == 0) || (it2.m > vrrmreq_max_.at(qam2)))
                        vrrmreq_max_[qam2] = it2.m;
                }
            }

            newmap[qam].push_back(vs);
            solvedquartets.insert(it);
        }

        PruneQuartets_(newtargets);

        targets = newtargets;
    } 

    // save the solved quartets
    for(const auto & it : solvedquartets)
        allam_.insert(it.amlist());
    allam_.insert({0,0,0,0});


    // remove all steps where the target m value is
    // the minimum.
    // This is called from within a loop over m
    // so it is handled by index math
    for(auto & it : newmap)
    {
        VRR_StepList vs = it.second;
        for(auto & it2 : vs)
        {
            if(it2.target.m == vrrmreq_min_[it.first])
                vrrmap_[it.first].push_back(it2);
        } 
    }


    // and empty step list for s s s s
    vrrmap_[{0,0,0,0}] = VRR_StepList();



    // determine requirements for all QAM
    // and the maximum integer constant
    maxint_ = 0;
    for(const auto & it : vrrmap_) // for all QAM
    {
        QAM qam = it.first;

        for(const auto & its : it.second)  // for all steps for this QAM
        {
            for(const auto & it3 : its.src)
            {
                if(it3)
                    qamreq_[qam].insert(it3.amlist());
            }

            std::string sstep = XYZStepToStr(its.xyz);

            maxint_ = std::max(its.ijkl[0], maxint_);
            maxint_ = std::max(its.ijkl[1], maxint_);
            maxint_ = std::max(its.ijkl[2], maxint_);
            maxint_ = std::max(its.ijkl[3], maxint_);

            if(its.type == RRStepType::I || its.type == RRStepType::J)
            {
                if(its.src[0])
                {
                    if(its.type == RRStepType::I)
                        varreq_[qam].insert(std::string("P_PA_") + sstep);
                    else
                        varreq_[qam].insert(std::string("P_PB_") + sstep);
                }

                if(its.src[1])
                    varreq_[qam].insert(std::string("aop_PQ_") + sstep);

                if(its.src[2] || its.src[3])
                {
                    varreq_[qam].insert("one_over_2p");
                    varreq_[qam].insert("a_over_p");
                    qamint_2p_[qam].insert(its.ijkl[0]);
                }

                if(its.src[4] || its.src[5])
                {
                    varreq_[qam].insert("one_over_2p");
                    varreq_[qam].insert("a_over_p");
                    qamint_2p_[qam].insert(its.ijkl[1]);
                }

                if(its.src[6])
                {
                    varreq_[qam].insert("one_over_2pq");
                    qamint_2pq_[qam].insert(its.ijkl[2]);
                }

                if(its.src[7])
                {
                    varreq_[qam].insert("one_over_2pq");
                    qamint_2pq_[qam].insert(its.ijkl[3]);
                }
            }
            else
            {
                if(its.src[0])
                {
                    if(its.type == RRStepType::K)
                        varreq_[qam].insert(std::string("Q_PA_") + sstep);
                    else
                        varreq_[qam].insert(std::string("Q_PB_") + sstep);

                }

                if(its.src[1])
                    varreq_[qam].insert(std::string("aoq_PQ_") + sstep);

                if(its.src[2] || its.src[3])
                {
                    varreq_[qam].insert("one_over_2q");
                    varreq_[qam].insert("a_over_q");
                    qamint_2q_[qam].insert(its.ijkl[2]);
                }

                if(its.src[4] || its.src[5])
                {
                    varreq_[qam].insert("one_over_2q");
                    varreq_[qam].insert("a_over_q");
                    qamint_2q_[qam].insert(its.ijkl[3]);
                }

                if(its.src[6])
                {
                    varreq_[qam].insert("one_over_2pq");
                    qamint_2pq_[qam].insert(its.ijkl[0]);
                }

                if(its.src[7])
                {
                    varreq_[qam].insert("one_over_2pq");
                    qamint_2pq_[qam].insert(its.ijkl[1]);
                }
            }
        }
    }


    // order the AM
    for(const auto & it : q) 
        AMOrder_AddWithDependencies_(amorder_, it.amlist());

}


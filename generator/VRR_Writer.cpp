#include "generator/ERIGeneratorInfo.hpp"
#include "generator/VRR_Writer.hpp"
#include "generator/Printing.hpp"
#include "generator/Naming.hpp"


///////////////////////////////
// Base VRR Writer
///////////////////////////////
VRR_Writer::VRR_Writer(const VRR_Algorithm_Base & vrr_algo)
    : vrr_algo_(vrr_algo)
{ 
}

bool VRR_Writer::HasVRR(void) const
{
    return vrr_algo_.HasVRR();
}

bool VRR_Writer::HasBraVRR(void) const
{
    return vrr_algo_.HasBraVRR();
}

bool VRR_Writer::HasKetVRR(void) const
{
    return vrr_algo_.HasKetVRR();
}

void VRR_Writer::DeclarePrimArrays(std::ostream & os, const ERIGeneratorInfo & info) const
{
    const VectorInfo & vinfo = info.GetVectorInfo();

    os << indent5 << "// Holds the auxiliary integrals ( i 0 | k 0 )^m in the primitive basis\n";
    os << indent5 << "// with m as the slowest index\n";


    for(const auto & qam : vrr_algo_.GetAllAM()) 
    {
        // add +1 fromm required m values to account for 0
        os << indent5 << "// AM = (" << qam[0] << " " << qam[1] << " | " << qam[2] << " " << qam[3] << " )\n";
        os << indent5 << vinfo.DoubleType() << " " << PrimVarName(qam)
           << "[" << (vrr_algo_.GetMReq(qam)+1) << " * " 
           << NCART(qam) << "] SIMINT_ALIGN_ARRAY_DBL;\n";
    }

    os << "\n\n";
}

void VRR_Writer::DeclarePrimPointers(std::ostream & os, const ERIGeneratorInfo & info) const
{
    for(const auto & qam : vrr_algo_.GetAllAM()) 
    {
        if(info.IsContQ(qam))
            os << indent4 << "double * restrict " << PrimPtrName(qam)
               << " = " << ArrVarName(qam) << " + abcd * " << NCART(qam) << ";\n";
    }

    os << "\n\n";
}

void VRR_Writer::WriteVRRSteps_(std::ostream & os, QAM qam, const VRR_StepSet & vs, const std::string & num_n,
                                const ERIGeneratorInfo & info) const
{
    const VectorInfo & vinfo = info.GetVectorInfo();

    os << indent5 << "// Forming " << PrimVarName(qam) << "[" << num_n << " * " << NCART(qam) << "];\n";

    os << indent5 << "for(n = 0; n < " << num_n << "; ++n)  // loop over orders of auxiliary function\n";
    os << indent5 << "{\n";

    os << "\n";

    // iterate over the requirements
    // should be in order since it's a set
    for(const auto & it : vs)
    {
        // Get the stepping
        XYZStep step = it.xyz;

        std::stringstream primname, srcname[8];
        primname << PrimVarName(qam) << "[n * " << NCART(qam) << " + " << it.target.idx() << "]";
        if(it.src[0])
            srcname[0] << PrimVarName(it.src[0].amlist()) << "[n * " << it.src[0].ncart() << " + " << it.src[0].idx() << "]";
        if(it.src[1])
            srcname[1] << PrimVarName(it.src[1].amlist()) << "[(n+1) * " << it.src[1].ncart() << " + " << it.src[1].idx() << "]";
        if(it.src[2])
            srcname[2] << PrimVarName(it.src[2].amlist()) << "[n * " << it.src[2].ncart() << " + " << it.src[2].idx() << "]";
        if(it.src[3])
            srcname[3] << PrimVarName(it.src[3].amlist()) << "[(n+1) * " << it.src[3].ncart() << " + " << it.src[3].idx() << "]";
        if(it.src[4])
            srcname[4] << PrimVarName(it.src[4].amlist()) << "[n * " << it.src[4].ncart() << " + " << it.src[4].idx() << "]";
        if(it.src[5])
            srcname[5] << PrimVarName(it.src[5].amlist()) << "[(n+1) * " << it.src[5].ncart() << " + " << it.src[5].idx() << "]";
        if(it.src[6])
            srcname[6] << PrimVarName(it.src[6].amlist()) << "[(n+1) * " << it.src[6].ncart() << " + " << it.src[6].idx() << "]";
        if(it.src[7])
            srcname[7] << PrimVarName(it.src[7].amlist()) << "[(n+1) * " << it.src[7].ncart() << " + " << it.src[7].idx() << "]";

        os << indent6 << "//" << it.target <<  " : STEP: " << step << "\n";

        std::string aoppq;
        std::string vrr_const0;
        std::string vrr_const1;
        std::string vrr_const2;
        std::string vrr_const3;
        std::string aover;

        // Note - the signs on a_over_p, a_over_q, etc, are taken care of in FileWriter.cpp
        if(it.type == RRStepType::I || it.type == RRStepType::J)
        {
            aoppq = std::string("aop_PQ_") + XYZStepToStr(step);
            aover = "a_over_p";
            vrr_const0 = StringBuilder("vrr_const_", it.ijkl[0], "_over_2p");
            vrr_const1 = StringBuilder("vrr_const_", it.ijkl[1], "_over_2p");
            vrr_const2 = StringBuilder("vrr_const_", it.ijkl[2], "_over_2pq");
            vrr_const3 = StringBuilder("vrr_const_", it.ijkl[3], "_over_2pq");
        }
        else
        {
            aoppq = std::string("aoq_PQ_") + XYZStepToStr(step);
            aover = "a_over_q";
            vrr_const0 = StringBuilder("vrr_const_", it.ijkl[2], "_over_2q");
            vrr_const1 = StringBuilder("vrr_const_", it.ijkl[3], "_over_2q");
            vrr_const2 = StringBuilder("vrr_const_", it.ijkl[0], "_over_2pq");
            vrr_const3 = StringBuilder("vrr_const_", it.ijkl[1], "_over_2pq");
        }



        if(info.HasFMA())
        {
            if(it.type == RRStepType::I || it.type == RRStepType::J)
            {
                os << indent6 << primname.str() << " = ";
                if(it.type == RRStepType::I)
                    os << "P_PA_" << step;
                else
                    os << "P_PB_" << step;

                os << " * " << srcname[0].str() << ";\n";
            }
            else
            {
                os << indent6 << primname.str();

                if(it.type == RRStepType::K)
                    os << " = Q_PA_" << step;
                else
                    os << " = Q_PB_" << step;

                os << " * " << srcname[0].str() << ";\n";
            }

            if(it.src[1])
                os << indent6 << primname.str() << " = " << vinfo.FMAdd(aoppq, srcname[1].str(), primname.str()) << ";\n"; 
            if(it.src[2] && it.src[3])
                os << indent6 << primname.str() << " = " << vinfo.FMAdd(vrr_const0, vinfo.FMAdd(aover, srcname[3].str(), srcname[2].str()), primname.str()) << ";\n"; 
            if(it.src[4] && it.src[5])
                os << indent6 << primname.str() << " = " << vinfo.FMAdd(vrr_const1, vinfo.FMAdd(aover, srcname[5].str(), srcname[4].str()), primname.str()) << ";\n"; 
            if(it.src[6])
                os << indent6 << primname.str() << " = " << vinfo.FMAdd(vrr_const2, srcname[6].str(), primname.str()) << ";\n";
            if(it.src[7])
                os << indent6 << primname.str() << " = " << vinfo.FMAdd(vrr_const3, srcname[7].str(), primname.str()) << ";\n";
        }
        else
        {
            if(it.type == RRStepType::I || it.type == RRStepType::J)
            {
                os << indent6 << primname.str();

                if(it.type == RRStepType::I)
                    os << " = P_PA_" << step;
                else
                    os << " = P_PB_" << step;

                os << " * " << srcname[0].str();
            }
            else
            {
                os << indent6 << primname.str();

                if(it.type == RRStepType::K)
                    os << " = Q_PA_" << step;
                else
                    os << " = Q_PB_" << step;

                os << " * " << srcname[0].str();
            }

            if(it.src[1])
                os << " + " << aoppq << " * " << srcname[1].str(); 
            if(it.src[2] && it.src[3])
                os << " + " << vrr_const0 << " * (" << srcname[2].str() << " + " << aover << " * " << srcname[3].str() << ")"; 
            if(it.src[4] && it.src[5])
                os << " + " << vrr_const1 << " * (" << srcname[4].str() << " + " << aover << " * " << srcname[5].str() << ")"; 
            if(it.src[6])
                os << " + " << vrr_const2 << " * " << srcname[6].str();
            if(it.src[7])
                os << " + " << vrr_const3 << " * " << srcname[7].str();
            os << ";\n";
        }
        
        os << "\n";

    }

    os << indent5 << "}\n";
}







///////////////////////////////////////
// Inline VRR Writer
///////////////////////////////////////

void VRR_Writer_Inline::AddConstants(ERIGeneratorInfo & info) const
{
    for(int i = 1; i <= vrr_algo_.GetMaxInt(); i++)
        info.AddIntegerConstant(i);
}


void VRR_Writer_Inline::WriteVRR(std::ostream & os, const ERIGeneratorInfo & info) const
{
    const VectorInfo & vinfo = info.GetVectorInfo();

    os << "\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << indent5 << "// Primitive integrals: Vertical recurrance\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << "\n";

    // constants
    for(const auto & it : vrr_algo_.GetAllInt_2p())
    {
        if(it == 1)
            os << indent5 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2p = one_over_2p;\n"; 
        else
            os << indent5 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2p = " << vinfo.IntConstant(it) << " * one_over_2p;\n"; 
    }

    for(const auto & it : vrr_algo_.GetAllInt_2q())
    {
        if(it == 1)
            os << indent5 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2q = one_over_2q;\n"; 
        else
            os << indent5 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2q = " << vinfo.IntConstant(it) << " * one_over_2q;\n"; 
    }

    for(const auto & it : vrr_algo_.GetAllInt_2pq())
    {
        if(it == 1)
            os << indent5 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2pq = one_over_2pq;\n"; 
        else
            os << indent5 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2pq = " << vinfo.IntConstant(it) << " * one_over_2pq;\n"; 
    }


    // iterate over increasing am
    for(const auto & qam : vrr_algo_.GetAMOrder())
    {
        // Write out the steps
        if(qam != QAM{0,0,0,0})
            WriteVRRSteps_(os, qam, vrr_algo_.GetSteps(qam), 
                           std::to_string(vrr_algo_.GetMReq(qam)+1),
                           info);

        os << "\n\n";
    }

    os << "\n\n";
}


void VRR_Writer_Inline::WriteVRRFile(std::ostream & os, std::ostream & osh, const ERIGeneratorInfo & info) const
{ }



///////////////////////////////////////
// External VRR Writer
///////////////////////////////////////

void VRR_Writer_External::AddConstants(ERIGeneratorInfo & info) const
{

}

void VRR_Writer_External::WriteVRR(std::ostream & os, const ERIGeneratorInfo & info) const
{
    os << "\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << indent5 << "// Primitive integrals: Vertical recurrance\n";
    os << indent5 << "//////////////////////////////////////////////\n";
    os << "\n";

    // iterate over increasing am
    for(const auto & am : vrr_algo_.GetAMOrder())
    {
        if(am != QAM{0,0,0,0})
        {
            os << indent5 << "VRR_" << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";
            os << indent7 << PrimVarName(am) << ",\n";

            for(const auto & it : vrr_algo_.GetAMReq(am))
                os << indent7 << PrimVarName(it) << ",\n";
            
            for(const auto & it : vrr_algo_.GetVarReq(am))
                os << indent7 << it << ",\n";

            os << indent7 << (vrr_algo_.GetMReq(am)+1) << ");\n";
            os << "\n";
        }
    }

    os << "\n\n";
}


void VRR_Writer_External::WriteVRRFile(std::ostream & os, std::ostream & osh, const ERIGeneratorInfo & info) const
{
    const VectorInfo & vinfo = info.GetVectorInfo();

    QAM am = info.FinalAM();

    os << "//////////////////////////////////////////////\n";
    os << "// VRR: ( " << amchar[am[0]] << " " << amchar[am[1]] << " | " << amchar[am[2]] << " " << amchar[am[3]] << " )\n";
    os << "//////////////////////////////////////////////\n";

    std::stringstream prototype;
    prototype << "void VRR_" << amchar[am[0]] << "_" << amchar[am[1]] << "_" << amchar[am[2]] << "_" << amchar[am[3]]  << "(\n";

    // final target
    prototype << indent3 << vinfo.DoubleType() << " * const restrict " << PrimVarName(am) << ",\n";

    for(const auto & it : vrr_algo_.GetAMReq(am))
        prototype << indent3 << vinfo.ConstDoubleType() << " * const restrict " << PrimVarName(it) << ",\n";
    
    for(const auto & it : vrr_algo_.GetVarReq(am))
        prototype << indent3 << vinfo.ConstDoubleType() << " " << it << ",\n";

    prototype << indent3 << "const int num_n)";

    os << prototype.str() << "\n";
    os << "{\n";


    // is this in standard order, or one of the permutations?
    if(am[1] == 0 && am[3] == 0)
    {
        // standard ( X s | Y s )
        // Make this one

        os << "    int n = 0;\n";

        for(const auto & it : vrr_algo_.GetIntReq_2p(am))
            os << indent1 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2p = " << vinfo.DoubleSet1(std::to_string(it)) << " * one_over_2p;\n"; 

        for(const auto & it : vrr_algo_.GetIntReq_2q(am))
            os << indent1 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2q = " << vinfo.DoubleSet1(std::to_string(it)) << " * one_over_2q;\n"; 

        for(const auto & it : vrr_algo_.GetIntReq_2pq(am))
            os << indent1 << vinfo.ConstDoubleType() << " vrr_const_" << it << "_over_2pq = " << vinfo.DoubleSet1(std::to_string(it)) << " * one_over_2pq;\n"; 

        // Write out the steps
        WriteVRRSteps_(os, am, vrr_algo_.GetSteps(am), "num_n", info);

    }
    else
    {
        // permutation. Call one of the others
        QAM tocall = am;
        if(am[0] == 0)
            std::swap(tocall[0], tocall[1]);
        if(am[2] == 0)
            std::swap(tocall[2], tocall[3]);

        os << indent1 << "// Routines are identical except for swapping of\n";
        os << indent1 << "// PA_x with PB_x and QC_x with QD_x\n";
        os << "\n";
        os << indent1 << "VRR_" << amchar[tocall[0]] << "_" << amchar[tocall[1]] << "_" << amchar[tocall[2]] << "_" << amchar[tocall[3]]  << "(\n";
    
        // final target
        os << indent3 << PrimVarName(am) << ",\n";
        
        for(const auto & it : vrr_algo_.GetAMReq(am))
            os << indent3 << PrimVarName(it) << ",\n";
    
        for(const auto & it : vrr_algo_.GetVarReq(am))
            os << indent3 << it << ",\n";

        os << indent3 << "num_n);";
    }

    os << "\n";
    os << "}\n";
    os << "\n\n";

    // header
    osh << prototype.str() << ";\n\n";
}

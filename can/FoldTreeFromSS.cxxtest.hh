// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>


#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

//package headers
#include <core/kinematics/Edge.hh>

/// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>


// C++ headers
#include <string>
#include <iostream>

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <utility/vector1.hh>



// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:


	void setUp() {
		core_init();
	}
    
    void test_hello_world(){
        TS_ASSERT("true");
    }
    
    utility::vector1< std::pair< core::Size, core::Size > >
    identify_secondary_structure_spans( std::string const & secstruct_codes ){
        utility::vector1< std::pair< core::Size, core::Size > > secondary_structure;
        core::Size start;
        for( core::Size iter = 0; iter < secstruct_codes.size()-1; ++iter){
            if (secstruct_codes[iter] =='H' || secstruct_codes[iter] == 'E'){
                start = iter;
                while (secstruct_codes[iter] == secstruct_codes[start]) {
                ++iter;
                }
        
            std::pair< core::Size, core::Size > newpair(start+1, iter);
            secondary_structure.push_back(newpair);
            --iter;
            }
        }//for
        return secondary_structure;
    }//end function

    void test_identify_secondary_structure_spans(){
        std::string str = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ";
        utility::vector1< std::pair< core::Size, core::Size > > test_result = identify_secondary_structure_spans(str);
        std::cout<<test_result<<std::endl;
        str = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE";
        test_result = identify_secondary_structure_spans(str);
        std::cout<<test_result<<std::endl;
        str="HHHHHHH   HHHHHHHHHHHH      HHHHHHHHHHHHEEEEEEEEEEHHHHHHH EEEEHHH ";
        test_result = identify_secondary_structure_spans(str);
        std::cout<<test_result;
        TS_ASSERT(true);
    }
    
    void fold_tree_from_ss(std::string secstruct_codes){
        //Initialize our foldtree
        core::kinematics::FoldTree foldtree;
        utility::vector1< std::pair< core::Size, core::Size > > test_result = identify_secondary_structure_spans(secstruct_codes);
        core::Size first_middle= (std::get<0>( test_result[1] ) +
                                  std::get<1> ( test_result[1]))/2;
        
        //Add the peptide edge 1st mid->start and 1st mid->1st end
        foldtree.add_edge(1, first_middle, core::kinematics::Edge::PEPTIDE);
        foldtree.add_edge(first_middle, std::get<1> (test_result[1]),
                          core::kinematics::Edge::PEPTIDE);
        
        //Add the all the jump edges and peptide edges 2-N-1
        for (core::Size iter = 2; iter < test_result.size(); ++iter) {
            core::Size cur_mid_sec = (std::get<0>( test_result[iter] ) +
                                  std::get<1> ( test_result[iter]))/2;
            core::Size cur_mid_loop = (std::get<1>( test_result[iter] ) +
                                       std::get<0> ( test_result[(iter+1)]))/2;
            
            //Add 1->iter jump edge
            foldtree.add_edge(first_middle, cur_mid_sec, 1);
            
            //determine if it is a real loop
            if ( cur_mid_loop != std::get<1>( test_result[iter] )
                && cur_mid_loop != std::get<0>( test_result[(iter+1)] ) ){
                add_edge( cur_mid_loop, std::get<1>( test_result[iter] ),
                         core::kinematics::Edge::PEPTIDE);
                add_edge ( cur_mid_loop, std::get<0>(test_result[(iter+1)] ) );
                
            }//endif
            
            core::Size last_mid = (std::get<0>( test_result[test_result.size()] ) +
                                   std::get<1> ( test_result[ test_result.size()]))/2;
            
            //add jump edge for 1->N
            add_edge( first_middle, last_mid, 1 );
            
            //add peptide edge for the last secondary structure
            add_edge( last_mid, std::get<0>( test_result[test_result.size()] ),
                     core::kinematics::Edge::PEPTIDE);
            add_edge( last_mid, secstruct_codes.size()+1,
                     core::kinematics::Edge::PEPTIDE);
            
                
        
        }
        
    }
	
};






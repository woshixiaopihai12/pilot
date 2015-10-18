//author=can_li
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>
#include <string>
#include <basic/options/option.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PymolMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/pack_rotamers.hh>





int
main(int argc, char ** argv) {
    try{
        devel::init( argc, argv );
        utility::vector1< std::string > filenames = basic::options::option[
                                                                           basic::options::OptionKeys::in::file::s ]();
        if ( filenames.size() > 0 ) {
            std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" <<
            std::endl;
            core::pose::PoseOP mypose = core::import_pose::pose_from_pdb( filenames[1] );
            core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();

            //intialize the MonteCarlo object
            protocols::moves::PyMolObserverOP observer = protocols::moves::AddPyMolObserver( *mypose, true, 0 );
            protocols::moves::MonteCarlo mc =  protocols::moves::MonteCarlo( *mypose, *sfxn, 0.5 );
            
            //set movemap
            core::kinematics::MoveMap mm;
            mm.set_bb( true );
            mm.set_chi( true );
            
            //set optimization
            core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
            
            //set atomtree minimizer
            core::optimization::AtomTreeMinimizer atm;
            
            //declare copy_pose
            core::pose::Pose copy_pose;

            for (core::Size iter = 1; iter<=1000; iter++) {
                std::cout<<"The score of last accepted pose :"<<mc.total_score_of_last_considered_pose()<<std::endl;
            //get the number of the residue to perturb
            core::Real uniform_random_number = numeric::random::uniform();
            core::Size randres = static_cast< core::Size >( mypose->total_residue() * uniform_random_number + 1 );
            //get the random number of perturbation
            core::Real pert1 = numeric::random::gaussian();
            core::Real pert2 = numeric::random::gaussian();
            //read the orginal phi and psi
            core::Real orig_phi = mypose->phi( randres );
            core::Real orig_psi = mypose->psi( randres );
            //perturb
            mypose->set_phi( randres, orig_phi + pert1 );
            mypose->set_psi( randres, orig_psi + pert2 );
            core::pack::task::PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( *mypose );
            repack_task->restrict_to_repacking();
            core::pack::pack_rotamers( *mypose, *sfxn, repack_task );
            copy_pose = *mypose;
            atm.run( copy_pose, mm, *sfxn, min_opts );
            *mypose = copy_pose;
            mc.boltzmann( *mypose );
                std::cout <<"Whether this mc move is accepted" <<mc.mc_accepted()<<std::endl;
            }
            
            
        } else {
            std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
            return 1;
        }
    } catch ( utility::excn::EXCN_Msg_Exception const & e ) {
        std::cout << "Caught exception! " << e.msg() << std::endl;
    }
    
    
} 

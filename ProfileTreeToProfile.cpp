/*---------------------------------------------------------------------------##
##  Library:
##      galosh::profuse
##  File:
##      ProfileTreeToProfile.cpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The profileTreeToProfile program.  It reads a ProfileTree from an .xml
##      file and converts _the root only_ to a Profile HMM (and saves it in a
##      profile file).
##
#******************************************************************************
#*
#*    This file is part of profuse, a suite of programs for working with
#*    Profile HMMs.  Please see the document CITING, which should have been
#*    included with this file.  You may use at will, subject to the license
#*    (Apache v2.0), but *please cite the relevant papers* in your documentation
#*    and publications associated with uses of this library.  Thank you!
#*
#*    Copyright (C) 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
#*    profuse is free software: you can redistribute it and/or modify it under
#*    the terms of the GNU Lesser Public License as published by the Free
#*    Software Foundation, either version 3 of the License, or (at your option)
#*    any later version.
#*
#*    profuse is distributed in the hope that it will be useful, but WITHOUT
#*    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#*    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for
#*    more details.
#*
#*    You should have received a copy of the GNU Lesser Public License along
#*    with profuse.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************/

#include "Algebra.hpp"
#include "ProfileTree.hpp"

#include <iostream>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

#include <seqan/basic.h>

#ifdef __HAVE_MUSCLE
int g_argc;
char **g_argv;
#endif // __HAVE_MUSCLE

using namespace galosh;

namespace galosh {

    template <class Serializable>
    static void
    readXML (
      Serializable & profuse_object,
      const char * filename
    )
    {
      // open the archive
      std::ifstream ifs( filename );
      assert( ifs.good() );
      boost::archive::xml_iarchive ia( ifs );
    
      // restore the profile from the archive
      ia >> BOOST_SERIALIZATION_NVP( profuse_object );
    } // readXML( Serializable &, const char * )

} // End namespace galosh

int
main ( int const argc, char const ** argv )
{
  typedef floatrealspace ProbabilityType;
#ifdef __PROFUSE_USE_AMINOS
  typedef seqan::AminoAcid20 ResidueType;
#else // __PROFUSE_USE_AMINOS .. else
  typedef seqan::Dna ResidueType;
#endif // __PROFUSE_USE_AMINOS .. else ..
  typedef ProfileTreeRoot<ResidueType, ProbabilityType> RootType;
  typedef ProfileTreeRoot<ResidueType, ProbabilityType> InternalNodeType;
  typedef ProfileTree<ResidueType, ProbabilityType, InternalNodeType > ProfileTreeType;

  if( argc < 2 ) {
    cout << "Usage: " << argv[ 0 ] << " <input ProfileTree.xml filename> [<output (galosh Profile) filename> ]" << endl;
    exit( 1 );
  }
  const bool be_verbose = true; //( argc >= 2 );
  if( be_verbose ) {
    cout << "Reading ProfileTree.xml from file '" << argv[ 1 ] << "'" << endl;
  }
  ProfileTreeType profile_tree;
  readXML( profile_tree, argv[ 1 ] );
  if( profile_tree.getProfileTreeRoot()->length() == 0 ) {
    cout << "No profiles were found in the ProfileTree.xml file '" << argv[ 1 ] << "'" << endl;
    return 1;
  } else if( profile_tree.nodeCount() > 1 ) {
    if( be_verbose ) {
      cout << "WARNING: Using only the root profile in the given ProfileTree.xml file." << endl;
    }
  }
  //if( be_verbose ) {
  //  cout << "\tgot:" << std::endl;
  //  cout << *( profile_tree.getProfileTreeRoot() );
  //  cout << endl;
  //}

  if( argc >= 3 ) {
    if( be_verbose ) {
      cout << "Writing Profile to file '" << argv[ 2 ] << "'" << endl;
    }
    std::ofstream profile_stream( argv[ 2 ] );
    assert( profile_stream.good() );
    profile_stream << *( profile_tree.getProfileTreeRoot() );
    profile_stream.close();
    if( be_verbose ) {
      cout << "\tdone." << endl;
    }
  } else {
    if( be_verbose ) {
      cout << "Profile is:" << endl;
    }
    cout << *( profile_tree.getProfileTreeRoot() );
    cout << endl;
  }

  exit( 0 );
} // main (..)

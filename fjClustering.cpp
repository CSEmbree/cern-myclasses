#include "fjClustering.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

fjClustering::fjClustering()
{
  std::cout<<" fjClustering called with default parameters "<<std::endl;
  double R = 1.0;
  *this = fjClustering(fastjet::kt_algorithm, R, fastjet::E_scheme, fastjet::Best);
}

fjClustering::fjClustering(fastjet::JetAlgorithm ja, 
			   double rparam,
			   fastjet::RecombinationScheme rs,
			   fastjet::Strategy s)
{
  rParameter = rparam;
  fjAlgorithm = ja;
  fjStrategy = s;
  fjRecombScheme = rs;
  fjJetDefinition = fastjet::JetDefinition(fjAlgorithm, rParameter, fjRecombScheme, fjStrategy);
  std::cout<<fjJetDefinition.description()<<std::endl;
}

fjClustering::~fjClustering()
{
  inputJets.clear();
  outputJets.clear();
 }

void fjClustering::ClearJets()
{
  outputJets.clear();
  //outputJets.resize(1, 0);
  inputJets.clear();
  //inputJets.resize(1, 0);
 }

void fjClustering::PrintJets()
{

  std::cout<<" input jets N="<<inputJets.size()<<std::endl;
  for (int  i = 0; i < inputJets.size(); i++)
  {
    std::cout<<i<<" px, py, py, E= "<<inputJets[i].px()
	                       <<" "<<inputJets[i].py()
	                       <<" "<<inputJets[i].pz()
	                       <<" "<<inputJets[i].e()
	     <<std::endl;
  }

  std::cout<<" output jets N= "<<outputJets.size()<<std::endl;
  for (int  i = 0; i < outputJets.size(); i++)
  {
    std::cout<<i<<" px, py, py, E= "<<outputJets[i].px()
                               <<" "<<outputJets[i].py()
                               <<" "<<outputJets[i].pz()
                               <<" "<<outputJets[i].e()
                               <<" id= "<<outputJets[i].user_index()
             <<std::endl;
  }


 }

void fjClustering::doClustering()
{

  fastjet::ClusterSequence cluster_seq(inputJets, fjJetDefinition);
  outputJets = cluster_seq.inclusive_jets(0.);
  //std::sort(cj.begin(), cj.end(), pT_sort());

}

void fjClustering::push_back(double px, double py, double pz, double E)
{
 //
  //fastjet::PseudoJet *jet= new fastjet::PseudoJet(px,py,pz,E);
  //inputJets.push_back(jet);
  inputJets.push_back(fastjet::PseudoJet(px,py,pz,E));

}

void fjClustering::push_back(double px, double py, double pz, double E, int pid)
{
 //
  //fastjet::PseudoJet *jet= new fastjet::PseudoJet(px,py,pz,E);
  //inputJets.push_back(jet);
  fastjet::PseudoJet *myjet2 = new fastjet::PseudoJet(px,py,pz,E);
  fastjet::PseudoJet myjet3 = *myjet2; 
  myjet3.set_user_index (pid);  
  inputJets.push_back(myjet3);
  //event.push_back(fastjet::PseudoJet(px,py,pz,E));
  delete myjet2;
  //  inputJets.push_back(fastjet::PseudoJet(px,py,pz,E));
  //inputJets.push_back(fastjet::PseudoJet(px,py,pz,E));

}


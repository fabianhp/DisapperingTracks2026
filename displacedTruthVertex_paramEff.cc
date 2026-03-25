//////////////////////////////////////////////////////
// ATLAS Displaced Vertex 13 TeV recast WITH TRUTH DV
// FOR LIGHT DISPLACED NEUTRALINO
// Written by Giovanna Cottin (gfcottin@gmail.com)
///////////////////////////////////////////////////
// ATLAS displaced vertex analysis in SUSY-2016-08
// Cuts: 
// Event Selection :
// Different triggers depending on analysis strategy
// DV Selection - Event pass if at least one displaced vertex with :
//   *   pT of tracks coming from DV  > 1 GeV
//   *   |d0| of tracks coming from DV  > 2 mm
//   This last two cuts define the signal region :
//   *   Number of tracks coming out of vertex > 5 
//   *   Invariant mass of vertex > 10 GeV 

#include "Pythia8/Pythia.h"
// Toy Detector Simulation (with input from Nishita Desai)
// Simulation includes tracklessjets by Giovanna Cottin
#include "ToyDetector-ATLAS-tracklessjet.h"
// NEW ATLAS parametrized_truthEff 
#include "parametrized_truthEff.h"


#include <vector>
#include <string>
#include <iterator>
#include <algorithm>

using namespace Pythia8;

int main(int argc, char *argv[]) {

  string a = argv[1]; // ENTER NEUTRALINO MASS FROM COMMAND LINE
  string b = argv[2]; // ENTER LAMBDA COUPLING FROM COMMAND LINE
  string c = argv[3];
  string d = argv[4];

  // You can always read an plain LHE file,
  // but if you ran "./configure --with-gzip" before "make"
  // then you can also read a gzipped LHE file.
  #ifdef GZIPSUPPORT
    bool useGzip = true;
  #else
    bool useGzip = false;
  #endif
    cout << " useGzip = " << useGzip << endl;

  Pythia pythia;
  Event& event = pythia.event;
  // Initialize Les Houches Event File run. List initialization information.
  pythia.readString("Beams:frameType = 4");
  /////////////////////////////////////////////////////////////////////////////////////

  pythia.readString("Beams:LHEF =/user/f/fhpinto/projects/MG5_aMC_v2_9_21/vquint-"+d+"/Events/run-"+a+"-"+b+"-"+c+"/unweighted_events.lhe.gz");
  



  int nEvent = 10000;//FROM LHE File
  int nAbort = 50;   

  ////// Initialize.
  pythia.init();

  //Files
  ofstream effCutFlow;
  effCutFlow.open("RPV_GRID/effCutFlow_GRID-"+a+"-"+b+"-"+c+".dat"); 
  ofstream decays;
  decays.open("RPV_GRID/decays_GRID-"+a+"-"+b+"-"+c+".dat"); 
  ofstream benchmarkParamsEff;
  benchmarkParamsEff.open("RPV_GRID/effParams_GRID-"+a+"-"+b+"-"+c+".dat"); 

  ofstream promptElec;
  promptElec.open("RPV_GRID/promptElec_GRID-"+a+"-"+b+"-"+c+".dat"); 
    
  // 3000 fb-1, 13 TeV
  // BR=1 

  //////////////////////////////////
  // width = 3.*1e11*6.581e-25/ctau
  //////////////////////////////////
  // Lifetime and mass checks
  cout<<"###################"<<endl;
  cout<<"pythia.particleData.tau0(9000005) = "<< setprecision(7)<< pythia.particleData.tau0(9000005)<<endl;
  cout<<"pythia.particleData.m0(9000006)   = "<< setprecision(7)<< pythia.particleData.m0(9000006)<<endl;
  cout<<"pythia.particleData.m0(9000007)   = "<< setprecision(7)<< pythia.particleData.m0(9000007)<<endl;


  int iAbort = 0;
  int nPromptElectron=0;
  int nDVReco=0;
  int nMaterial=0;
  int nFidutial=0;
  int nTrk=0;
  int nMass=0;
  int nDVEff=0;

  ToyDetector detector;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; ; ++iEvent) {
    cout << "######## Event "<<iEvent<<" #########"<<endl;
  // for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }
    if(!detector.getObjects(event)) {
        cout << "No objects found" << endl;
      } 
      else {
         // cout<<"Size of muons = "<< detector.muons.size()<<endl;
         // cout << "Size of jets      = "<< detector.jets.size() <<endl;
          cout << "Size of electrons = "<< detector.electrons.size() <<endl;
         // detector.printObjects();
      }
    /////////////////
    //Event Selection
    /////////////////
    // Cuts on prompt electrons
      /////////////////////////////////////////////
      //Displaced Vertex truth identification
      //Find all particles coming from a NEUTRALINO
      /////////////////////////////////////////////  
      std::vector<int> motherIndices;   
//      //Truth neutralino mass test
//      for (int i= 0; i < event.size(); i++){
//        if (abs(event[i].id()) == 1000022){
//          double mass = event[i].m();
//          cout<<"MASS = "<<mass <<endl;
//        }
//      }
      //Get neutralinos from event
      for (int i= 0; i < event.size(); i++){
        if ((abs(event[i].id()) == 9000005) || (abs(event[i].id()) == 9000006)){
          cout<<"There is a v+/v++"<<endl;
          double trux = event[i].xDec() - event[i].xProd();
          double truy = event[i].yDec() - event[i].yProd();
          double truz = event[i].zDec() - event[i].zProd();    
          double vProd = sqrt(pow2(event[i].xProd())+pow2(event[i].yProd())+pow2(event[i].zProd()));     
           //If neutralino is displaced 
          if ((abs(trux) > 0.00001) || (abs(truy) > 0.00001)) {
            double decayLength = sqrt(pow2(trux) + pow2(truy) + pow2(truz));
            double p = sqrt(pow2(event[i].px())+pow2(event[i].py())+pow2(event[i].pz()));
            double mass = event[i].m();
            cout<<"Neutralino MASS = "<<mass <<endl;
            double gamma =  p/mass;
            // Pythia stores energy units in GeV, with c=1. But for tau0, is in units of mm/c
            double dist = event[i].vDec().pAbs();//should be equal to tau0*betagamma and to decayLenght
            cout<<"###################"<<endl;
            cout<<"dist (from vDec) [mm] =  "<< event[i].vDec().pAbs()<<endl;
            cout<<"decayLength [mm]         =  "<< decayLength<<endl;
            cout<<"gamma                =  "<< gamma<<endl;

            //In this case dcyLen and dist are the same. vDec() gives the decay vertex coordinates,  
            //calculated from the production vertex (at origin) and the proper lifetime 
            //cout << "decay vertex coordinates are = "<<event[i].vDec()<<endl;
            // decays<<pythia.particleData.tau0(6000012) <<" "<<gamma<<" "<<decayLength<<endl;
            //decays<<pythia.particleData.tau0(1000022) <<" "<<gamma<<" "<<decayLength<<endl;
          }
        }    
    }
  }      //End of event loop.
    ///////////////////////////////////////////////////////////////
    // Cross section (needs to be multiplied by overall efficiency)
    ///////////////////////////////////////////////////////////////
    double width_p;
    double width_pp;
    if((pythia.particleData.tau0(9000005)>0)){ 
        width_pp = 3.*1e11*6.581e-25/pythia.particleData.tau0(9000005);
        width_p = 3.*1e11*6.581e-25/pythia.particleData.tau0(9000006);      
    }
    double sigma = pythia.info.sigmaGen(); // in mb
    cout << "####################################################################" << endl;
    cout << " Cross section is "  << sigma * 1e12 << " [fb] "<<endl;
    cout << "####################################################################" << endl;
    ////////////////////////
    // CutFlow default ATLAS
    ////////////////////////
    // Number of events - Relative Efficiency - Overall efficiency
   // effCutFlow <<" All events "        << nEvent                << " " << nEvent*100./nEvent            << " " << nEvent*100./nEvent    << endl;
   // effCutFlow <<" Prompt Elec  "      << nPromptElectron       << " " << nPromptElectron*100./nEvent   << " " << nPromptElectron*100./nEvent<< endl;
   // effCutFlow <<" DV Fidutial "       << nFidutial             << " " << nFidutial*100./nPromptElectron<< " " << nFidutial*100./nEvent << endl;
   // effCutFlow <<" DV Ntrk "           << nTrk                  << " " << nTrk*100./nFidutial           << " " << nTrk*100./nEvent      << endl;
   // effCutFlow <<" DV Mass "           << nMass                 << " " << nMass*100./nTrk               << " " << nMass*100./nEvent     << endl;
   // effCutFlow <<" DV Eff "            << nDVEff                << " " << nDVEff*100./nMass             << " " << nDVEff*100./nEvent    << endl;
    
    benchmarkParamsEff<<setprecision(9)<<fixed<<pythia.particleData.tau0(9000005)<<" "<<pythia.particleData.tau0(9000006)<<" "<<width_p<<" "<<width_pp<<setprecision(9)<<fixed<<" "<<nDVEff<<" "<<nDVEff*1.0/nEvent<<" "<<a<<" "<<b<<" "<<c<<" "<<sigma*1e12<<endl;

    ////////////////////////    
   // cout <<" CutFlow for ctau, tau = "<< pythia.particleData.tau0(1000022)<<" , "<<pythia.particleData.tau0(1000022)/(299.792)    <<endl;
   // cout <<" All events "        << nEvent                << " " << nEvent*100./nEvent            << " " << nEvent*100./nEvent    << endl;
  //  cout <<" Prompt Elec  "      << nPromptElectron       << " " << nPromptElectron*100./nEvent   << " " << nPromptElectron*100./nEvent<< endl;
  //  cout <<" DV Fidutial "       << nFidutial             << " " << nFidutial*100./nPromptElectron<< " " << nFidutial*100./nEvent << endl;
   // cout <<" DV Ntrk "           << nTrk                  << " " << nTrk*100./nFidutial           << " " << nTrk*100./nEvent      << endl;
    //cout <<" DV Mass "           << nMass                 << " " << nMass*100./nTrk               << " " << nMass*100./nEvent     << endl;
   //cout <<" DV Eff "            << nDVEff                << " " << nDVEff*100./nMass             << " " << nDVEff*100./nEvent    << endl;
  
    pythia.stat();
    effCutFlow.close();  
    benchmarkParamsEff.close();
    decays.close();
    promptElec.close();
    //massNtrk.close();
  return 0;
}


#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"

int materialCheck() {
    TFile* file = TFile::Open("geant4_material_tracks.root");
    TTree* tree = (TTree*)file->Get("material-tracks");

    tree->Print();
    tree->Draw("mat_x:mat_y:mat_z", "mat_z<16567", "");
    return 0;
}
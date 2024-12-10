////////////////////////////////////////////////////////////////////////
// Global parameters
////////////////////////////////////////////////////////////////////////

TString geoTag = "Prototype";
TString geoPath;
TString GTPCGeoParamsFile;
Double_t fHalfSizeTPC_X;
Double_t fHalfSizeTPC_Y;
Double_t fHalfSizeTPC_Z;
Double_t fSizeOfVirtualPad;
Double_t fMaxDriftTime;
Int_t histoBins;
Int_t histoBins2;
char hname[255];

void initializeGlobals() {
    geoPath = gSystem->Getenv("VMCWORKDIR");
    GTPCGeoParamsFile = geoPath + "/glad-tpc/params/HYDRAprototype_FileSetup.par";
    GTPCGeoParamsFile.ReplaceAll("//", "/");

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BGTPCGeoPar* geoPar = (R3BGTPCGeoPar*)rtdb->getContainer("GTPCGeoPar");
    if (!geoPar) {
        cout << "No R3BGTPCGeoPar can be loaded from the rtdb" << endl;
        exit(1);
    }
    R3BGTPCGasPar* gasPar = (R3BGTPCGasPar*)rtdb->getContainer("GTPCGasPar");
    if (!gasPar) {
        cout << "No R3BGTPCGasPar can be loaded from the rtdb" << endl;
        exit(1);
    }

    FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo(); // Ascii file
    parIo1->open(GTPCGeoParamsFile, "in");
    rtdb->setFirstInput(parIo1);
    rtdb->initContainers(0);

    fHalfSizeTPC_X = geoPar->GetActiveRegionx() / 2.; // 50cm in X (row)
    fHalfSizeTPC_Y = geoPar->GetActiveRegiony() / 2.; // 20cm in Y (time)
    fHalfSizeTPC_Z = geoPar->GetActiveRegionz() / 2.; // 100cm in Z (column)
    fSizeOfVirtualPad = geoPar->GetPadSize();         // 1: pads of 1cm^2 , 10: pads of 1mm^2
    fMaxDriftTime = round((geoPar->GetActiveRegiony() / gasPar->GetDriftVelocity()) * pow(10, -3)); // us

    histoBins = 2 * fHalfSizeTPC_X * fSizeOfVirtualPad;
    histoBins2 = 2 * fHalfSizeTPC_Z * fSizeOfVirtualPad;
}

void laserPadPlane(TString inputSimFile, TString title, TString mode)
{
    initializeGlobals();

    TH2D* htrackInPads = 0;
    TH2D* hdriftTimeInPads = 0;
    TH2D* hdepth1InPads = 0;
    TH2D* hdepth2InPads = 0;
    TH1S** h1_ProjPoint_TimeExample = 0;

    htrackInPads = new TH2D("htrackInPads_"+title,
                            "All tracks in the XZ Pads Plane",
                            histoBins, 0, histoBins, histoBins2, 0, histoBins2); // in [pad number]
    htrackInPads->SetYTitle("Z [pad number]");
    htrackInPads->SetXTitle("X [pad number]");



    hdriftTimeInPads = new TH2D("hdriftTimeInPads"+title,
                                "All tracks in the XZ Pads Plane with drift time",
                                histoBins,
                                0,
                                histoBins,
                                histoBins2,
                                0,
                                histoBins2); // in [pad number]
    hdriftTimeInPads->SetYTitle("Z [pad number]");
    hdriftTimeInPads->SetXTitle("X [pad number]");

    hdepth1InPads = new TH2D("hdepth1InPads"+title,
                             "track In the Drift-Z Pads Plane",
                             histoBins,
                             0,
                             fMaxDriftTime,
                             histoBins2,
                             0,
                             histoBins2);
    hdepth1InPads->SetYTitle("Z [pad number]");
    hdepth1InPads->SetXTitle("(drift) time [us]");

    hdepth2InPads = new TH2D("hdepth2InPads"+title,
                             "track In the Drift-X Pads Plane",
                             histoBins,
                             0,
                             fMaxDriftTime,
                             histoBins,
                             0,
                             histoBins);
    hdepth2InPads->SetYTitle("X [pad number]");
    hdepth2InPads->SetXTitle("(drift) time [us]");

    ////////////////////////////////////////////////////////////////////////
    // EVENT ACCESS
    ////////////////////////////////////////////////////////////////////////

    // File access
    TFile* simFile = TFile::Open(inputSimFile, "READ");
    TTree* TEvt = (TTree*)simFile->Get("evt");
    Int_t nevents = TEvt->GetEntries();

    // Projection Point Definition
    TClonesArray* gtpcProjPointCA;
    R3BGTPCProjPoint* ppoint = new R3BGTPCProjPoint;
    gtpcProjPointCA = new TClonesArray("R3BGTPCProjPoint", 5);
    TBranch* branchGTPCProjPoint = TEvt->GetBranch("GTPCProjPoint");
    branchGTPCProjPoint->SetAddress(&gtpcProjPointCA);

    // Event 0 analysis
    Int_t padsPerEvent = 0;
    Int_t nb = 0;
    Int_t beamPadsWithSignalPerEvent, productPadsWithSignalPerEvent;
    Double_t tPad;
    Int_t xPad, zPad, yPad;
    Int_t numberOfTimeHistos = 0;

    gtpcProjPointCA->Clear();
    nb += TEvt->GetEvent(0);
    padsPerEvent = gtpcProjPointCA->GetEntries();

    if (padsPerEvent > 0)
    {
        h1_ProjPoint_TimeExample = new TH1S*[padsPerEvent];
        for (Int_t h = 0; h < padsPerEvent; h++)
        {
            ppoint = (R3BGTPCProjPoint*)gtpcProjPointCA->At(h);

            tPad = ((TH1S*)(ppoint->GetTimeDistribution()))->GetMean();
            htrackInPads->GetBinXYZ(ppoint->GetVirtualPadID(), xPad, zPad, yPad);
            xPad--;
            zPad--;

            htrackInPads->Fill(xPad, zPad, ppoint->GetCharge());
            hdriftTimeInPads->Fill(xPad, zPad, tPad);
            hdepth1InPads->Fill(tPad, zPad, ppoint->GetCharge());
            hdepth2InPads->Fill(tPad, xPad, ppoint->GetCharge());

            sprintf(hname, "pad %i", ppoint->GetVirtualPadID());
            h1_ProjPoint_TimeExample[h] = (TH1S*)((ppoint->GetTimeDistribution()))->Clone(hname);
        }

        numberOfTimeHistos = padsPerEvent;
    }

    ////////////////////////////////////////////////////////////////////////
    // RESULTS
    ////////////////////////////////////////////////////////////////////////

    auto *file = new TFile("laser_results.root", mode);
    gROOT->SetStyle("Default");
    htrackInPads->Draw("zcol");
    htrackInPads->Write();

    file->Close();


}


void laserPadPlane2(TString inputSimFile, TString title, TString mode)
{
    initializeGlobals();

    std::shared_ptr<R3BGTPCMap> fTPCMap;
    TH2* fPadPlane;
    fTPCMap = std::make_shared<R3BGTPCMap>();
    fTPCMap->GeneratePadPlane();
    fPadPlane = fTPCMap->GetPadPlane();


    ////////////////////////////////////////////////////////////////////////
    // EVENT ACCESS
    ////////////////////////////////////////////////////////////////////////

    // File access
    TFile* simFile = TFile::Open(inputSimFile, "READ");
    TTree* TEvt = (TTree*)simFile->Get("evt");
    Int_t nevents = TEvt->GetEntries();

    // Projection Point Definition
    TClonesArray* gtpcProjPointCA;
    R3BGTPCProjPoint* ppoint = new R3BGTPCProjPoint;
    gtpcProjPointCA = new TClonesArray("R3BGTPCProjPoint", 5);
    TBranch* branchGTPCProjPoint = TEvt->GetBranch("GTPCProjPoint");
    branchGTPCProjPoint->SetAddress(&gtpcProjPointCA);

    // Event 0 analysis
    Int_t padsPerEvent = 0;
    Int_t nb = 0;
    Int_t beamPadsWithSignalPerEvent, productPadsWithSignalPerEvent;
    Double_t tPad;
    Int_t xPad, zPad, yPad;
    Int_t numberOfTimeHistos = 0;

    gtpcProjPointCA->Clear();
    nb += TEvt->GetEvent(0);
    padsPerEvent = gtpcProjPointCA->GetEntries();

    if (padsPerEvent > 0)
    {
        for (Int_t h = 0; h < padsPerEvent; h++)
        {
            ppoint = (R3BGTPCProjPoint*)gtpcProjPointCA->At(h);

            tPad = ((TH1S*)(ppoint->GetTimeDistribution()))->GetMean();
            if (ppoint->GetVirtualPadID() < 0){continue;}
            std::cout << ppoint->GetVirtualPadID() << std::endl;
            
            fPadPlane->GetBinXYZ(ppoint->GetVirtualPadID(), xPad, zPad, yPad);
            fPadPlane->Fill(xPad, zPad);
            //xPad--;
            //zPad--;
            fPadPlane->AddBinContent(xPad, zPad, ppoint->GetCharge());

            sprintf(hname, "pad %i", ppoint->GetVirtualPadID());
        }

        numberOfTimeHistos = padsPerEvent;
    }

    ////////////////////////////////////////////////////////////////////////
    // RESULTS
    ////////////////////////////////////////////////////////////////////////

    auto *file = new TFile("laser_results.root", mode);
    fPadPlane->Draw("zcol");
    fPadPlane->Write();

    file->Close();


}

/*
void fillMagneticField(TH3F* &hx, TH3F* &hy, TH3F* &hz, Double_t fYIn=24.312,
 Double_t fZIn=6.98, Double_t fAlpha=11.67, Double_t fBeta=1.36)
{

    initializeGlobals();

    // Some tpc parameters
    Double_t fHalfSizeTPC_X = 4.4;            
    Double_t fHalfSizeTPC_Y = 14.7;     
    Double_t fHalfSizeTPC_Z = 12.8; 
    Double_t fXIn = 0.;
    Double_t TargetOffsetX =  40;
    Double_t TargetOffsetY = 0;
    Double_t TargetOffsetZ = 260;      
    Double_t TargetOffsetZ_FM = 263.4; 
    Double_t TargetAngle = 14. * TMath::Pi() / 180;

    ////////////////////////////////////////////////////////////////////////
    // MAGNETIC FIELD INITIALIZATION
    ////////////////////////////////////////////////////////////////////////

    R3BGladFieldMap *gladField = new R3BGladFieldMap("R3BGladMap");
    gladField->SetScale(-1.);
    gladField->Init();

    ////////////////////////////////////////////////////////////////////////
    // LASER LIMITS
    ////////////////////////////////////////////////////////////////////////

    // Express the angles in rad 
	fAlpha *= TMath::Pi() / 180.;
	fBeta *= TMath::Pi() / 180.;

    // Calculate the radius where the laser scapes the TPC
    Double_t rX = 2 * fHalfSizeTPC_X / cos(fBeta) / sin(fAlpha);
    Double_t rY, rZ;
    
    if (fAlpha>0) // The laser only can scape through the top
    {
        rY = (2 * fHalfSizeTPC_Y - fYIn) / sin(fBeta);
    }
    else // The laser only can scape through the bottom
    {
        rY = (fYIn) / sin(fBeta);
    }

    
    if ((fBeta < TMath::Pi() / 2) && (fBeta > 0.))
    {
        rZ = (2 * fHalfSizeTPC_Z - fZIn) / cos(fBeta) / cos(fAlpha);
    }
    else
    {
        rZ = (fZIn) / cos(fBeta) / cos(fAlpha);
    }

    Double_t rads[3] = {rX, rY, rZ}; 
    Double_t rMin = TMath::MinElement(3, rads);

    // Get the range
    Double_t fXOut = rMin * cos(fBeta) * sin(fAlpha);
    Double_t fYOut = fZIn + rMin * cos(fBeta) * cos(fAlpha);
    Double_t fZOut = fYIn + rMin * sin(fBeta);

    Double_t ele_x_init, ele_y_init, ele_z_init;
    Double_t B_x, B_y, B_z;
    Double_t padX, padZ;
    TH3F *h = (TH3F*)hx->Clone("h");
    Int_t histoBins = hx->GetNbinsX();
    Int_t histoBins2 = hx->GetNbinsY();
    Double_t xPad, zPad, yPad;
    
    for (Double_t xval = 0; xval < 2*fHalfSizeTPC_X; xval+=.5)
    {
        for (Double_t yval = 0; yval < 2*fHalfSizeTPC_Y; yval+=.5)
        {
            for (Double_t zval = 0; zval < 2*fHalfSizeTPC_X; zval+=.5)
            {
                ele_y_init = yval -fHalfSizeTPC_Y;
                ele_x_init =+ cos(-TargetAngle) * (xval) + sin(-TargetAngle) * (zval);
                ele_z_init = (TargetOffsetZ_FM - fHalfSizeTPC_Z) - sin(-TargetAngle) * (xval) + cos(-TargetAngle) * (zval);

                B_x = 0.1 * gladField->GetBx(ele_x_init, ele_y_init, ele_z_init); 
                B_y = 0.1 * gladField->GetBy(ele_x_init, ele_y_init, ele_z_init);
                B_z = 0.1 * gladField->GetBz(ele_x_init, ele_y_init, ele_z_init);

                padX = histoBins * xval / 2. / fHalfSizeTPC_X;
		        padZ = histoBins2 * zval / 2. / fHalfSizeTPC_Z;


                Int_t hBin = h->Fill(padX, padZ, yval);
                h->GetBinXYZ(hBin, xPad, zPad, yPad);
                xPad--;
                yPad--;
                zPad--;

                std::cout << xPad << " " << zPad << " " << yPad << " " << B_x << std::endl;
                std::cout << xPad << " " << zPad << " " << yPad << " " << B_y << std::endl;
                std::cout << xPad << " " << zPad << " " << yPad << " " << B_z << std::endl;
                std::cout << " " << std::endl;


                hx->SetBinContent(xPad, zPad, yPad, B_x);
                hy->SetBinContent(xPad, zPad, yPad, B_y);
                hz->SetBinContent(xPad, zPad, yPad, B_z);

            }
        }
    }



}
*/



void laserVis(TString gladFieldSource = "../proj/laser_gen_gladField.root", TString constFieldSource = "no_init", Bool_t plotMagField=kFALSE)
{

    initializeGlobals();

    // Save the pad plane projections using the former function
    laserPadPlane2(gladFieldSource, "gladField", "RECREATE");
    if (constFieldSource != "no_init"){laserPadPlane2(constFieldSource, "constField", "UPDATE");}
    
    // Get both histograms
    auto *f = new TFile("laser_results.root", "read");
    auto *histo1 = (TH2S*)f->Get("htrackInPads_constField"); 
    auto *histo2 = (TH2S*)f->Get("htrackInPads_gladField"); 

    
    //if (plotMagField)
    //{
        /*
        auto *hx = new TH3F("hx", "X Coordinate", histoBins, 0, histoBins,
                                             histoBins2, 0, histoBins2,
                                             100, 0, 100);

        */
        /*
        auto *hx = new TH3F("hx", "X Coordinate", 500, 0, 44,
                                             500, 0, 120,
                                             200, 0, 200);
        

        auto *hy = (TH3F*)hx->Clone("hy");
        auto *hz = (TH3F*)hx->Clone("hz");

        fillMagneticField(hx, hy, hz);

        hz->Draw("zcol");
        hy->Draw("zcol");
        hx->Draw("zcol");
 

    }
    */


    if (histo1 == NULL)
    {
        histo2->Draw();


        ////// TERMINAR DE PONER ESTO BIEN :))))
    }
    else
    {
        // Get the projections
        auto *projConst = histo1->ProfileX();
        auto *projGlad = histo2->ProfileX();

        // Draw and compare both results
        auto *c = new TCanvas();
        projConst->Draw();
        projGlad->SetLineColor(kRed);
        projGlad->Draw("same");
    }




}
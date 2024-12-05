void run_laser_gen(Bool_t constB = kFALSE)
{

    TStopwatch timer;
    timer.Start();

    // Input file: simulation
    TString inFile = "../sim/Prototype/sim.root";
    // Input file: parameters
    TString parFile = "../sim/Prototype/par.root";
    // Output file
    TString outFile = "laser_gen";


    // -----   Create analysis run   ----------------------------------------
    FairRunAna* fRun = new FairRunAna();

    // ----------------- Magnetic Field Implementation -------------------------
    FairField *magField = NULL;

    if (constB)
    {
        magField = new R3BFieldConst();
        double B_y = 100.; //[kG]
        ((R3BFieldConst*) magField)->SetField(0., B_y, 0.);
        ((R3BFieldConst*) magField)->SetFieldRegion(-200.0, // x_min
                                200.0,  // x_max
                                -100.0, // y_min
                                100.0,  // y_max
                                -150.0, // z_min
                                450.0); // z_max

        outFile += "_constField.root";

    }
    else
    {
        magField = new R3BGladFieldMap("R3BGladMap");
        ((R3BGladFieldMap*) (magField))->SetScale(-1.);
        ((R3BGladFieldMap*) (magField))->Init();

        outFile += "_gladField.root";
        
    }

    fRun->SetField(magField);
    fRun->SetSource(new FairFileSource(inFile));
    fRun->SetOutputFile(outFile.Data());

    // -----   Runtime database   ---------------------------------------------
    FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
    FairParRootFileIo* parIn = new FairParRootFileIo(kTRUE);
    parIn->open(parFile.Data());
    rtdb->setFirstInput(parIn);
    rtdb->print();

    // -------------- Laser Generation -----------------------------------------
    R3BGTPCLaserGen* lasergen = new R3BGTPCLaserGen();
    lasergen->SetDriftParameters(15.e-9, 0.0048, 0.00000216, 0.00000216, 2);
    //lasergen->SetLaserParameters(9.71, 1.54, 0., 16.204, 5.76, 20000., 0.1);
    lasergen->SetLaserParameters(11.67, 1.36, 0., 24.312, 6.98, 2);
    //lasergen->SetOutputMode(1);
    lasergen->SetSizeOfVirtualPad(5); // 1 means pads of 1cm^2, 10 means pads of 1mm^2, ...
    lasergen->SetNumberOfGeneratedElectrons(10000);

    fRun->AddTask(lasergen);
    fRun->Init();
    fRun->Run(0, 0);
    delete fRun;

    timer.Stop();

    cout << "Macro finished succesfully!" << endl;
    cout << "Output file writen: " << outFile << endl;
    cout << "Parameter file writen: " << parFile << endl;
    cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << endl;
}

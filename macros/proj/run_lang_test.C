void run_lang_test()
{
    TStopwatch timer;
    timer.Start();

    // Input file: simulation
    TString inFile = "../sim/Prototype/sim.root";
    // Input file: parameters
    TString parFile = "../sim/Prototype/par.root";
    // Output file
    TString outFile = "lang_test.root";




    // -----   Create analysis run   ----------------------------------------
    FairRunAna* fRun = new FairRunAna();

    /*
    R3BFieldConst *constField = new R3BFieldConst();

    double B_y = 21.; //[kG]
    (constField)->SetField(0., B_y, 0.);
    (constField)->SetFieldRegion(-200.0, // x_min
                             200.0,  // x_max
                             -100.0, // y_min
                             100.0,  // y_max
                             -150.0, // z_min
                             450.0); // z_max

    fRun->SetField(constField);
    */

    R3BGladFieldMap *magField = new R3BGladFieldMap("R3BGladMap");
    magField->SetScale(-1.);
    magField->Init();
    fRun->SetField(magField);



    fRun->SetSource(new FairFileSource(inFile));
    fRun->SetOutputFile(outFile.Data());

    // -----   Runtime database   ---------------------------------------------
    FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
    FairParRootFileIo* parIn = new FairParRootFileIo(kTRUE);
    parIn->open(parFile.Data());
    rtdb->setFirstInput(parIn);
    rtdb->print();

    R3BGTPCLangevinTest* lantest = new R3BGTPCLangevinTest();
    lantest->SetDriftParameters(15.e-9, 0.0048, 0.00000216, 0.00000216, 2);
    //lantest->SetLaserParameters(9.71, 1.54, 0., 16.204, 5.76);
    lantest->SetLaserParameters(11.67, 1.36, 0., 24.312, 6.98, 5000);
    lantest->SetSizeOfVirtualPad(5); // 1 means pads of 1cm^2, 10 means pads of 1mm^2, ...
    lantest->SetNumberOfGeneratedElectrons(10000);

    fRun->AddTask(lantest);

    fRun->Init();
    fRun->Run(0, 0);
    delete fRun;

    timer.Stop();

    cout << "Macro finished succesfully!" << endl;
    cout << "Output file writen: " << outFile << endl;
    cout << "Parameter file writen: " << parFile << endl;
    cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << endl;
}

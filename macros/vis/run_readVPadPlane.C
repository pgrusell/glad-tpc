{
    gROOT->ProcessLine(".L readVPadPlane.C");
    reader("../proj/lang_test.root"); // change the folder according to the detector you want to use
    guiForPads(0);
}

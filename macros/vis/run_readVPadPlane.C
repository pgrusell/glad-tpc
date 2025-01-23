{
    gROOT->ProcessLine(".L readVPadPlane.C");
    reader("../proj/Prototype/lang.root"); // change the folder according to the
                                           // detector you want to use
    guiForPads(0);
}

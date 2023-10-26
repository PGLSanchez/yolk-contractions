setAutoThreshold("Default no-reset");
for (i = 1; i <= nSlices; i++) {
    setSlice(i);
    run("Create Selection");
    run("Measure");
}

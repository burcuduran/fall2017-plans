To run the script we need three input argument
First Load the script by writing ( .L heepcheck_plot.C)
Then do heepcheck("target","sign","hpart");
There are three target. So you can put target1 or target2 or target3
For sign you can put add or subs or no, to add/subtract or ignore energyloss calcuation.
For hpart (HMS particle) you can put e or p.
You also need the input file heepcheck.input which has the following information:
RunNumber, Beam Energy, e-angle, e-momentum, p-angle, p-momentum, beam_offset, e-angle offset, e-momentum offset, p-angle offset, p-momentum offset
Only the last line of the input file is used, all other lines are ignored.
We have included an example input file and root file from E94139, run number 23460, which can be used to run and test the script.



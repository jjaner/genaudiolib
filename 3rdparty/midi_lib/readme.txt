Midi Tools
http://www.ee.columbia.edu/~csmit/matlab_midi.html


Extending the midi toolbox in matlab

The Midi Toolbox is a great tool for analyzing midi files, but it can't handle files with tempo changes properly. I've also had trouble using it with Windows 7. My code can substitute for the read_midi.m and write_midi.m functions of the Midi Toolbox. 
Installation:

Download the midi jar KaraokeMidiJava.jar.
Open classpath.txt (type "edit classpath.txt" in Matlab) and add a line containing the location of KaraokeMidiJava.jar. Save the file* and restart Matlab.
Download and unzip the matlab code midi_lib.zip. Add the location of this directory to the Matlab path (File -> Set Path...)*.
Examples:

Simple reading and writing: example_script1.html.
Extensions: example_script2.html.

* If you don't have write access to classpath.txt or pathdef.m (the file that is modified by File -> Set Path...), you can make local copies of these files, which you can modify. Then, when you start Matlab, you will need to specify the location of these files in the command line. For example, in Linux/Unix, you would use the command "Matlab /path/to/classpath.txt /path/to/pathdef.m."


package cellfatereprogramming;

import fileOperations.DeleteFile;
import java.util.ArrayList;
import javax.script.ScriptException;
import stablemotifs.*;

/**
 *
 * @author Jorge G. T. Zañudo

Copyright (c) 2015 Jorge G. T. Zañudo and Réka Albert.
 
The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */ 

public class TLGLNetwork {
    
        public static void main(String[] args) throws ScriptException {
        
        String fileName="TLGLNetwork.txt";    
        long startTime,estimatedTime;
        String[] simplifiedFunctions;
        String[] names;
        String networkName;
        String sinksFilename,transitionsFilename;
        ArrayList<String[]> finalAttractorsReduction;
        ArrayList<ArrayList<String[]>> resultNoDuplicates;
        ArrayList<String[]> unstableAttractorsFinal;
        ArrayList<String[]> stableMotifs;
        ArrayList<String[][]> attractors;
        ArrayList<String[]> attractorsResult;
        ArrayList<String[]> motifsAttractorsDiagram;
        ArrayList<String[]> controlSets;
        ArrayList<String> controlSetAttractors;
        String[] motifsToRemove;
        String[][] attractorsToMerge;
        
        Network nw;
        SuccessionDiagram succeDiag;
        
        networkName=fileName.split("\\.")[0];
        System.out.println("\nFilename: "+fileName);
        System.out.println("Creating Boolean table directory: "+networkName);
        ReadWriteFiles.createTablesFromBooleanRules(networkName, fileName);
        System.out.println("Boolean table directory created.");
        System.out.println("Creating functions and names files.");
        nw=OtherMethods.RecreateNetwork(networkName);
        nw.findNodeOutputs(); 
        names=nw.getNames();
        simplifiedFunctions=nw.getFunctions();
        System.out.println("Functions and names files created.");
        startTime = System.nanoTime();
        System.out.println("Performing network reduction...");
        attractors=NetworkReduction.FullNetworkReductionTop(names, simplifiedFunctions,networkName);
        System.out.println("Network reduction complete.");
        estimatedTime = System.nanoTime() - startTime;       
        System.out.println("Removing duplicate quasi-attractors.");
        attractorsResult=new ArrayList<String[]>();        
        for(int i=0;i<attractors.size();i++){attractorsResult.add(attractors.get(i)[1]);}
        resultNoDuplicates=NetworkReduction.removeDuplicateQuasiattractors(attractorsResult,names);
        finalAttractorsReduction=resultNoDuplicates.get(0);
        unstableAttractorsFinal=resultNoDuplicates.get(1);
        System.out.println("Total number of quasi-attractors: "+(finalAttractorsReduction.size()+unstableAttractorsFinal.size()));
        System.out.println("Number of putative quasi-attractors: "+unstableAttractorsFinal.size());
        System.out.println("Total time for finding quasi-attractors: "+((double) estimatedTime)/((double) 1000000000)+" s");
        System.out.println("Writing TXT files with quasi-attractors and stable motifs.");
        stableMotifs=ReadWriteFiles.getStableMotifsFromFileAndRemoveDuplicates("StableMotifs-"+networkName+".txt");
        ReadWriteFiles.writeStableMotifsAndQuasiAttractorFiles(networkName,names,stableMotifs,finalAttractorsReduction,unstableAttractorsFinal);
        DeleteFile.deletefile("StableMotifs-"+networkName+".txt");
        NetworkReduction.writeReducedNetwork(finalAttractorsReduction, networkName, "QA", names, simplifiedFunctions);
        NetworkReduction.writeReducedNetwork(unstableAttractorsFinal, networkName, "PQA", names, simplifiedFunctions);
        System.out.println("Starting analyis of stable motif succession diagram.");
        startTime = System.nanoTime();
        System.out.println("Identifying quasi-attractors corresponding to stable motif sequences.");
        motifsAttractorsDiagram=NetworkReduction.getAttractorsCorrespondingToMotifSequence(finalAttractorsReduction, names, networkName);
        ReadWriteFiles.writeAttractorsCorrespondingToMotifSequence(motifsAttractorsDiagram,networkName);
        DeleteFile.deletefile("DiagramSinks-"+networkName+".txt");
        sinksFilename=networkName+"-DiagramSequencesAttractors.txt";
        transitionsFilename="Diagram-"+networkName+".txt";
        if(networkName.equals("TLGLNetwork")){
            motifsToRemove=new String[2];motifsToRemove[0]="(P2=0)";motifsToRemove[1]="(P2=1)";
            attractorsToMerge=new String[2][3];
            attractorsToMerge[0][0]="Attractor2";attractorsToMerge[1][0]="Leukemia";
            attractorsToMerge[0][1]="Attractor1";attractorsToMerge[1][1]="Leukemia";
            attractorsToMerge[0][2]="Attractor0";attractorsToMerge[1][2]="Apoptosis";            
            OtherMethods.simplifySequences(sinksFilename, transitionsFilename, motifsToRemove, attractorsToMerge);
            transitionsFilename="Diagram-"+networkName+"Modified.txt";
            sinksFilename=networkName+"-DiagramSequencesAttractorsModified.txt";
        }
        if(networkName.equals("ThNetwork")){
            motifsToRemove=new String[0];
            attractorsToMerge=new String[2][12];
            attractorsToMerge[0][0]="Attractor0";attractorsToMerge[1][0]="Th17";
            attractorsToMerge[0][1]="Attractor10";attractorsToMerge[1][1]="Th2";
            attractorsToMerge[0][2]="Attractor11";attractorsToMerge[1][2]="Th2";
            for(int i=1;i<=3;i++){attractorsToMerge[0][i+2]="Attractor"+i;attractorsToMerge[1][i+2]="Treg";}
            for(int i=4;i<=9;i++){attractorsToMerge[0][i+2]="Attractor"+i;attractorsToMerge[1][i+2]="Th1";}            
            OtherMethods.simplifySequences(sinksFilename, transitionsFilename, motifsToRemove, attractorsToMerge);
            transitionsFilename="Diagram-"+networkName+"Modified.txt";
            sinksFilename=networkName+"-DiagramSequencesAttractorsModified.txt";
        }
        succeDiag=new SuccessionDiagram (sinksFilename,transitionsFilename,names,simplifiedFunctions);
        succeDiag.findControlSets();
        estimatedTime = System.nanoTime() - startTime;
        System.out.println("Total time for finding stable motif control sets: "+((double) estimatedTime)/((double) 1000000000)+" s");
        System.out.println("Writing TXT files with stable motif control sets.");
        succeDiag.writeStableMotifControlSets(networkName+"-StableMotifControlSets.txt");
        System.out.println("Done!");
        controlSets=succeDiag.getControlSets();
        controlSetAttractors=succeDiag.getControlSetAttractors();
        System.out.println("\nCONTROL SETS: ");
        for(int i=0;i<controlSets.size();i++){
            for(String str: controlSets.get(i)){
                System.out.print(str+" ");
            }
            System.out.print("-> "+controlSetAttractors.get(i)+"\n");
        }
        

    }

    
}

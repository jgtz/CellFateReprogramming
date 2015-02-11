package cellfatereprogramming;

import fileOperations.FileToRead;
import java.util.*;
import javax.script.ScriptException;

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


public class MinimalSubsetsOfStableMotifsThNetwork {
    
    public static void main(String[] args) throws ScriptException{
        
        String fileName="ThNetwork-DiagramSequencesAttractorsModified.txt";
        FileToRead fr=new FileToRead(fileName);
        String line,attractor;
        String[] lineArray,motifs;
        ArrayList<String> attractors=new ArrayList<>();
        ArrayList<String> attractorsCurrent,attractorsNext,attractorsFinal;
        Map<String,Integer> map=new HashMap();
        Map<Integer,String> mapReverse=new HashMap();
        ArrayList<HashSet> setArray=new ArrayList<>();
        ArrayList<HashSet> setCurrent,setNext,setFinal;
        ArrayList<Integer> list;
        HashSet set;
        boolean bool,boolSave,boolSubset;
        Iterator itr;


        map.put("(GATA3=0-IL4=0-IL4R_2=0-STAT6=0)", 0);
        map.put("(IFNG=0-IFNGR=0-STAT1=0-TBET=0)", 1);
        map.put("(IL23=1-IL23R=1-STAT3=1)", 2);
        map.put("(FOXP3=1)", 3);
        map.put("(FOXP3=0)", 4);
        map.put("(IL10=1-IL10R=1-STAT3=1)", 5);
        map.put("(IL21=1-IL21R=1-STAT3=1)", 6);
        map.put("(RORGT=1)", 7);
        map.put("(RORGT=0)", 8);
        map.put("(IL10=0-IL10R=0-IL21=0-IL21R=0-IL23R=0-STAT3=0)", 9);
        map.put("(TBET=1)", 10);
        map.put("(TBET=0)", 11);
        map.put("(GATA3=0)", 12);
        map.put("(GATA3=1)", 13);
        map.put("(GATA3=0-TBET=1)", 14);
        map.put("(GATA3=1-TBET=0)", 15);
        map.put("(FOXP3=0-IFNG=1-IFNGR=1-STAT1=1)", 16);
        mapReverse.put(0,"(GATA3=0-IL4=0-IL4R_2=0-STAT6=0)");
        mapReverse.put(1,"(IFNG=0-IFNGR=0-STAT1=0-TBET=0)");
        mapReverse.put(2,"(IL23=1-IL23R=1-STAT3=1)");
        mapReverse.put(3,"(FOXP3=1)");
        mapReverse.put(4,"(FOXP3=0)");
        mapReverse.put(5,"(IL10=1-IL10R=1-STAT3=1)");
        mapReverse.put(6,"(IL21=1-IL21R=1-STAT3=1)");
        mapReverse.put(7,"(RORGT=1)");
        mapReverse.put(8,"(RORGT=0)");
        mapReverse.put(9,"(IL10=0-IL10R=0-IL21=0-IL21R=0-IL23R=0-STAT3=0)");
        mapReverse.put(10,"(TBET=1)");
        mapReverse.put(11,"(TBET=0)");
        mapReverse.put(12,"(GATA3=0)");
        mapReverse.put(13,"(GATA3=1)");
        mapReverse.put(14,"(GATA3=0-TBET=1)");
        mapReverse.put(15,"(GATA3=1-TBET=0)");
        mapReverse.put(16,"(FOXP3=0-IFNG=1-IFNGR=1-STAT1=1)");

        for(int i=0;i<17;i++){
            System.out.println("MOTIF "+i+": "+mapReverse.get((Integer) i));
        }


        while(fr.hasNext()){
            line=fr.nextLine();
            lineArray=line.split("\t");
            motifs=lineArray[0].split(" ");
            set=new HashSet<>();
            bool=true;
            for(int i=0;i<motifs.length;i++){
                set.add(map.get(motifs[i]));
            }
            for(int i=0;i<setArray.size();i++){
                if(setArray.get(i).equals(set)){
                    bool=false;
                    break;
                }
            }
            if(bool){
                setArray.add(set);
                attractors.add(lineArray[1]);
            }

        }

        fr.close();
        setFinal=new ArrayList<>();
        attractorsFinal=new ArrayList<>();
        setCurrent=new ArrayList<>();
        attractorsCurrent=new ArrayList<>();
        for(int i=0;i<setArray.size();i++){
            setCurrent.add(setArray.get(i));
            attractorsCurrent.add(attractors.get(i));
        }

        while(!attractorsCurrent.isEmpty()){
            setNext=new ArrayList<>();
            attractorsNext=new ArrayList<>();
            for(int i=0;i<setCurrent.size();i++){
                attractor=attractorsCurrent.get(i);
                list=new ArrayList<>(setCurrent.get(i));
                boolSave=true;
                for(int j=0;j<list.size();j++){
                    boolSubset=true;
                    set=new HashSet();
                    for(int k=0;k<list.size();k++){
                        if(k!=j){set.add(list.get(k));}
                    }
                    for(int k=0;k<setArray.size();k++){
                        if(setArray.get(k).containsAll(set)){
                            if(!attractors.get(k).equals(attractor)){
                                boolSubset=false;
                                break;
                            }
                        }
                    }
                    if(boolSubset){
                        boolSave=false;
                        setNext.add(set);
                        attractorsNext.add(attractor);
                    }
                }
                if(boolSave){
                    bool=true;
                    for(int j=0;j<setFinal.size();j++){
                        if(setFinal.get(j).equals(setCurrent.get(i))){
                            bool=false;
                            break;
                        }
                    }
                    if(bool){
                        setFinal.add(setCurrent.get(i));
                        attractorsFinal.add(attractor);
                    }
                }
            }
            setCurrent=new ArrayList<>();
            attractorsCurrent=new ArrayList<>();
            for(int i=0;i<setNext.size();i++){
                setCurrent.add(setNext.get(i));
                attractorsCurrent.add(attractorsNext.get(i));
            }
        }
        
        System.out.println("\nMinimal subsets of stable motifs that are sufficient for a sequence to lead to a single helper T cell subtype");
        System.out.print("Th1: ");
        for(int i=0;i<setFinal.size();i++){        
            if(attractorsFinal.get(i).equals("Th1")){
                System.out.print(setFinal.get(i).toString()+" ");
            }
        }
        System.out.print("\nTh2: ");
        for(int i=0;i<setFinal.size();i++){        
            if(attractorsFinal.get(i).equals("Th2")){
                System.out.print(setFinal.get(i).toString()+" ");
            }
        }
        System.out.print("\nTh17: ");
        for(int i=0;i<setFinal.size();i++){     
            if(attractorsFinal.get(i).equals("Th17")){
                System.out.print(setFinal.get(i).toString()+" ");
            }
        }
        System.out.print("\nTreg: ");
        for(int i=0;i<setFinal.size()-1;i++){    
            if(attractorsFinal.get(i).equals("Treg")){           
                System.out.print(setFinal.get(i).toString()+" ");
            }
        }
        System.out.print("\n");
 
    }
    
}

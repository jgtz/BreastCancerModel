/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bcnetwork;

import booleandynamicmodeling.Network;
import booleandynamicmodeling.OtherMethods;
import booleandynamicmodeling.ReadWriteFiles;
import booleandynamicmodeling.UpdateMethods;
import fileOperations.FileToWrite;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import javax.script.ScriptException;

/**
 *
 * @author Jorge G. T. Zañudo
 * @date November 2017
 */
public class NetworkDoublePerturbationSimulations {

    //This program simulates the ER+ breast cancer network model in the presence of each possible
    //double node perturbation, starting from the cancer steady state initial condition. 
    
    //The regulatory functions of the modelbare given in the TXT file "BreastCancerModel.txt".
    //The model dynamics are governed by the stochastic general asynchronous updating
    //scheme. Every time step corresponds to several (Ntimes) updates, with one time
    //step corresponding to the average number of updates needed to update a slow node. 

     /**
     * @param args args[0] is the name of the TXT file where the model is. For the breast cancer model it is "BreastCancerModel.txt".
     */
    
        public static void main(String[] args) throws ScriptException {
        
        String fileName="BreastCancerModel.txt"; //This file contains the Boolean rules of the model       
        ArrayList<Integer> fast=new ArrayList<>(); //This array stores the indices of the fast nodes
        ArrayList<Integer> slow=new ArrayList<>(); //This array stores the indices of the slow nodes
        double Apofraction1,Apofraction2,Apofraction3,Prolfraction1,Prolfraction2,Prolfraction3,Prolfraction4,Prolfraction,Apofraction;
        int Apofraction1Ind,Apofraction2Ind,Apofraction3Ind,Prolfraction1Ind,Prolfraction2Ind,Prolfraction3Ind,Prolfraction4Ind;
        int N;
        int Ntime; //One time step is Ntime updates 
        int IC=10000; //Number of initial conditions
        int T=100; //This is the number of time steps
        
        double p; //This is the probability of updating any of the slow nodes at any update step
        //the update probability of a fast node is 5 times that of a slow node
        int[] nodeStates,pastState;
        int updateNode,ranState;       
        int index;
        String line="";
        int input,inputState,input2,inputState2;
        ArrayList<Integer> state1=new ArrayList<Integer>();state1.add(1);state1.add(0);
        ArrayList<Integer> state2=new ArrayList<Integer>();state2.add(1);state2.add(0);
        FileToWrite fw=new FileToWrite("BreastCancerDoublePerturbations.txt");
        //The average apoptosis and proliferation at the end of each simulation for each
        //double perturbation will be saved in this tab-separated file
                
        System.out.println("\nFilename: "+fileName);
        System.out.println("Creating Boolean table directory: "+fileName.split("\\.")[0]);
        ReadWriteFiles.createTablesFromBooleanRules(fileName.split("\\.")[0], fileName);
        System.out.println("Boolean table directory created.");
        System.out.println("Creating functions and names files.");
        Network nw=OtherMethods.RecreateNetwork(fileName.split("\\.")[0]);        
        System.out.println("Functions and names files created.");
        nw.findNodeOutputs(); 
        N=nw.getN();
        HashMap namesDictionary=new HashMap<Integer,String>();
        HashMap indexDictionary=new HashMap<String,Integer>();
        for(int i=0;i<N;i++){namesDictionary.put(i, nw.getNames()[i]);indexDictionary.put(nw.getNames()[i],i);}
        nodeStates=new int[N];
        Apofraction1Ind=(int) indexDictionary.get("Apoptosis");Apofraction2Ind=(int) indexDictionary.get("Apoptosis_2");Apofraction3Ind=(int) indexDictionary.get("Apoptosis_3");
        Prolfraction1Ind=(int) indexDictionary.get("Proliferation");Prolfraction2Ind=(int) indexDictionary.get("Proliferation_2");Prolfraction3Ind=(int) indexDictionary.get("Proliferation_3");Prolfraction4Ind=(int) indexDictionary.get("Proliferation_4");
        
        
        //Slow nodes for the simulations. The nodes not added to the "slow" list are the fast nodes and will stay in the "fast" list
        for(int i=0;i<N;i++){fast.add(new Integer(i));}
        fast.remove(new Integer((int)indexDictionary.get("HER3")));slow.add(new Integer((int)indexDictionary.get("HER3")));
        fast.remove(new Integer((int)indexDictionary.get("HER3_2")));slow.add(new Integer((int)indexDictionary.get("HER3_2")));
        fast.remove(new Integer((int)indexDictionary.get("BIM")));slow.add(new Integer((int)indexDictionary.get("BIM")));
        fast.remove(new Integer((int)indexDictionary.get("BCL2")));slow.add(new Integer((int)indexDictionary.get("BCL2")));
        fast.remove(new Integer((int)indexDictionary.get("MCL1")));slow.add(new Integer((int)indexDictionary.get("MCL1")));
        fast.remove(new Integer((int)indexDictionary.get("HER2")));slow.add(new Integer((int)indexDictionary.get("HER2")));
        fast.remove(new Integer((int)indexDictionary.get("HER3_T")));slow.add(new Integer((int)indexDictionary.get("HER3_T")));
        fast.remove(new Integer((int)indexDictionary.get("cycE_CDK2_T")));slow.add(new Integer((int)indexDictionary.get("cycE_CDK2_T")));
        fast.remove(new Integer((int)indexDictionary.get("IGF1R_T")));slow.add(new Integer((int)indexDictionary.get("IGF1R_T")));
        fast.remove(new Integer((int)indexDictionary.get("BCL2_T")));slow.add(new Integer((int)indexDictionary.get("BCL2_T")));
        fast.remove(new Integer((int)indexDictionary.get("BIM_T")));slow.add(new Integer((int)indexDictionary.get("BIM_T")));
        fast.remove(new Integer((int)indexDictionary.get("ER")));slow.add(new Integer((int)indexDictionary.get("ER")));
        fast.remove(new Integer((int)indexDictionary.get("ESR1")));slow.add(new Integer((int)indexDictionary.get("ESR1")));
        fast.remove(new Integer((int)indexDictionary.get("ESR1_2")));slow.add(new Integer((int)indexDictionary.get("ESR1_2")));
        fast.remove(new Integer((int)indexDictionary.get("FOXA1")));slow.add(new Integer((int)indexDictionary.get("FOXA1")));
        fast.remove(new Integer((int)indexDictionary.get("PBX1")));slow.add(new Integer((int)indexDictionary.get("PBX1")));
        fast.remove(new Integer((int)indexDictionary.get("ER_transcription")));slow.add(new Integer((int)indexDictionary.get("ER_transcription")));
        fast.remove(new Integer((int)indexDictionary.get("ER_transcription_2")));slow.add(new Integer((int)indexDictionary.get("ER_transcription_2")));
        fast.remove(new Integer((int)indexDictionary.get("cyclinD")));slow.add(new Integer((int)indexDictionary.get("cyclinD")));
        fast.remove(new Integer((int)indexDictionary.get("cyclinD_2")));slow.add(new Integer((int)indexDictionary.get("cyclinD_2")));
        fast.remove(new Integer((int)indexDictionary.get("CDK46")));slow.add(new Integer((int)indexDictionary.get("CDK46")));
        fast.remove(new Integer((int)indexDictionary.get("E2F")));slow.add(new Integer((int)indexDictionary.get("E2F")));
        fast.remove(new Integer((int)indexDictionary.get("E2F_2")));slow.add(new Integer((int)indexDictionary.get("E2F_2")));
        fast.remove(new Integer((int)indexDictionary.get("E2F_3")));slow.add(new Integer((int)indexDictionary.get("E2F_3")));
        fast.remove(new Integer((int)indexDictionary.get("p21_p27_T")));slow.add(new Integer((int)indexDictionary.get("p21_p27_T")));
        fast.remove(new Integer((int)indexDictionary.get("SGK1_T")));slow.add(new Integer((int)indexDictionary.get("SGK1_T")));
        fast.remove(new Integer((int)indexDictionary.get("PDK1")));slow.add(new Integer((int)indexDictionary.get("PDK1")));
        fast.remove(new Integer((int)indexDictionary.get("ER")));slow.add(new Integer((int)indexDictionary.get("ER")));
        fast.remove(new Integer((int)indexDictionary.get("mTORC2")));slow.add(new Integer((int)indexDictionary.get("mTORC2")));
        fast.remove(new Integer((int)indexDictionary.get("FOXO3_Ub")));slow.add(new Integer((int)indexDictionary.get("FOXO3_Ub")));
        fast.remove(new Integer((int)indexDictionary.get("IGF1R")));slow.add(new Integer((int)indexDictionary.get("IGF1R")));
        fast.remove(new Integer((int)indexDictionary.get("IGF1R_2")));slow.add(new Integer((int)indexDictionary.get("IGF1R_2")));
        fast.remove(new Integer((int) indexDictionary.get("Apoptosis")));slow.add(new Integer((int) indexDictionary.get("Apoptosis")));
        fast.remove(new Integer((int) indexDictionary.get("Apoptosis_2")));slow.add(new Integer((int) indexDictionary.get("Apoptosis_2")));
        fast.remove(new Integer((int) indexDictionary.get("Apoptosis_3")));slow.add(new Integer((int) indexDictionary.get("Apoptosis_3")));
        fast.remove(new Integer((int) indexDictionary.get("Proliferation")));slow.add(new Integer((int) indexDictionary.get("Proliferation")));
        fast.remove(new Integer((int) indexDictionary.get("Proliferation_2")));slow.add(new Integer((int) indexDictionary.get("Proliferation_2")));
        fast.remove(new Integer((int) indexDictionary.get("Proliferation_3")));slow.add(new Integer((int) indexDictionary.get("Proliferation_3")));
        fast.remove(new Integer((int) indexDictionary.get("Proliferation_4")));slow.add(new Integer((int) indexDictionary.get("Proliferation_4")));
        fast.remove(new Integer((int)indexDictionary.get("MYC")));slow.add(new Integer((int)indexDictionary.get("MYC")));
        fast.remove(new Integer((int)indexDictionary.get("MYC_2")));slow.add(new Integer((int)indexDictionary.get("MYC_2")));
        
        p=1.0*slow.size()/(1.0*slow.size()+5.0*fast.size());//p=pslow, 5*p=pfast, p*slow.size()+pfast*fast.size()= 1 =p(slow.size()+5*fast.size())      
        Ntime=slow.size()+5*fast.size();
        Ntime=(int)(Ntime/5.0);
        System.out.println("Perturbation\tApofrac1\tApofrac2\tApofrac3\tApofrac\tProlfrac1\tProlfrac2\tProlfrac3\tProlfrac4\tProlfrac");                
        fw.writeLine("Perturbation\tApofrac1\tApofrac2\tApofrac3\tApofrac\tProlfrac1\tProlfrac2\tProlfrac3\tProlfrac4\tProlfrac");
        
        for(Integer s1: state1){
            for(Integer s2: state2){
                for(int i1=0;i1<N;i1++){
                    for(int i2=0;i2<i1;i2++){
                        input=i1;
                        inputState=s1.intValue();
                        input2=i2;
                        inputState2=s2.intValue();
                        Apofraction1=0;Apofraction2=0;Apofraction3=0;Prolfraction1=0;Prolfraction2=0;Prolfraction3=0;Prolfraction4=0;Apofraction=0;Prolfraction=0;

                        for(int r=0;r<IC;r++){
                            for(int i=0;i<N;i++){nodeStates[i]=0;}
                                //Nodes state of sources
                                nodeStates[(int) indexDictionary.get("IGF1R_T")]=1;
                                nodeStates[(int) indexDictionary.get("HER2")]=0;
                                nodeStates[(int) indexDictionary.get("HER3_T")]=0;
                                nodeStates[(int) indexDictionary.get("PDK1")]=0;
                                nodeStates[(int) indexDictionary.get("SGK1_T")]=0;
                                nodeStates[(int) indexDictionary.get("mTORC2")]=1;
                                nodeStates[(int) indexDictionary.get("PIM")]=0;
                                nodeStates[(int) indexDictionary.get("PTEN")]=0;

                                nodeStates[(int) indexDictionary.get("Fulvestrant")]=0;
                                nodeStates[(int) indexDictionary.get("Alpelisib")]=1;
                                nodeStates[(int) indexDictionary.get("Everolimus")]=0;
                                nodeStates[(int) indexDictionary.get("Palbociclib")]=0;
                                nodeStates[(int) indexDictionary.get("Trametinib")]=0;

                                ranState=(int)(1.5*Math.random());
                                nodeStates[(int) indexDictionary.get("BIM")]=ranState;
                                nodeStates[(int) indexDictionary.get("BIM_T")]=ranState;
                                nodeStates[(int) indexDictionary.get("BAD")]=0;
                                nodeStates[(int) indexDictionary.get("MCL1")]=1;           
                                nodeStates[(int) indexDictionary.get("BCL2_T")]=ranState;
                                nodeStates[(int) indexDictionary.get("BCL2")]=ranState;

                                //Nodes state of cancer attractor            
                                nodeStates[(int) indexDictionary.get("IGF1R")]=1;
                                nodeStates[(int) indexDictionary.get("IGF1R_2")]=0;
                                nodeStates[(int) indexDictionary.get("HER2_3")]=0;
                                nodeStates[(int) indexDictionary.get("SGK1")]=0;
                                nodeStates[(int) indexDictionary.get("RAS")]=1;
                                nodeStates[(int) indexDictionary.get("RAS_2")]=0;
                                nodeStates[(int) indexDictionary.get("MAPK")]=1;
                                nodeStates[(int) indexDictionary.get("MAPK_2")]=0;
                                nodeStates[(int) indexDictionary.get("PI3K")]=1;
                                nodeStates[(int) indexDictionary.get("PIP3")]=1;
                                nodeStates[(int) indexDictionary.get("PDK1_pm")]=1;
                                nodeStates[(int) indexDictionary.get("mTORC2_pm")]=1;
                                nodeStates[(int) indexDictionary.get("AKT")]=1;
                                nodeStates[(int) indexDictionary.get("p21_p27_T")]=0;
                                nodeStates[(int) indexDictionary.get("p21_p27")]=0;
                                nodeStates[(int) indexDictionary.get("cycE_CDK2")]=1;
                                nodeStates[(int) indexDictionary.get("cycE_CDK2_T")]=1;            
                                nodeStates[(int) indexDictionary.get("KMT2D")]=0;
                                nodeStates[(int) indexDictionary.get("TSC")]=0;
                                nodeStates[(int) indexDictionary.get("PRAS40")]=0;
                                nodeStates[(int) indexDictionary.get("mTORC1")]=1;
                                nodeStates[(int) indexDictionary.get("FOXO3")]=0;
                                nodeStates[(int) indexDictionary.get("FOXO3_Ub")]=0;
                                nodeStates[(int) indexDictionary.get("EIF4F")]=1;
                                nodeStates[(int) indexDictionary.get("S6K")]=1;
                                nodeStates[(int) indexDictionary.get("Translation")]=1;
                                nodeStates[(int) indexDictionary.get("ER")]=1;
                                nodeStates[(int) indexDictionary.get("ESR1")]=1;
                                nodeStates[(int) indexDictionary.get("ESR1_2")]=0;
                                nodeStates[(int) indexDictionary.get("FOXA1")]=0;
                                nodeStates[(int) indexDictionary.get("PBX1")]=1;
                                nodeStates[(int) indexDictionary.get("ER_transcription")]=1;
                                nodeStates[(int) indexDictionary.get("ER_transcription_2")]=0;
                                nodeStates[(int) indexDictionary.get("cyclinD")]=1;
                                nodeStates[(int) indexDictionary.get("cyclinD_2")]=0;
                                nodeStates[(int) indexDictionary.get("CDK46")]=1;
                                nodeStates[(int) indexDictionary.get("cycD_CDK46")]=1;
                                nodeStates[(int) indexDictionary.get("cycD_CDK46_2")]=0;
                                nodeStates[(int) indexDictionary.get("pRb")]=1;
                                nodeStates[(int) indexDictionary.get("pRb_2")]=1;
                                nodeStates[(int) indexDictionary.get("pRb_3")]=0;
                                nodeStates[(int) indexDictionary.get("E2F")]=1;
                                nodeStates[(int) indexDictionary.get("E2F_2")]=1;
                                nodeStates[(int) indexDictionary.get("E2F_3")]=0;
                                nodeStates[(int) indexDictionary.get("Proliferation")]=1;
                                nodeStates[(int) indexDictionary.get("Proliferation_2")]=1;
                                nodeStates[(int) indexDictionary.get("Proliferation_3")]=1;
                                nodeStates[(int) indexDictionary.get("Proliferation_4")]=0;
                                nodeStates[(int) indexDictionary.get("Apoptosis")]=0;
                                nodeStates[(int) indexDictionary.get("Apoptosis_2")]=0;
                                nodeStates[(int) indexDictionary.get("Apoptosis_3")]=0;
                                nodeStates[(int) indexDictionary.get("MYC")]=1;
                                nodeStates[(int) indexDictionary.get("MYC_2")]=0;


                            for(int t=0;t<T;t++){
                                if(t==0){nodeStates[input]=inputState;nodeStates[input2]=inputState2;}
                                for(int n=0;n<Ntime;n++){
                                    if(Math.random()>=p){index=(int)(Math.random()*fast.size());updateNode=fast.get(index);}
                                    else{index=(int)(Math.random()*slow.size());updateNode=slow.get(index);}
                                    pastState=Arrays.copyOf(nodeStates, nodeStates.length);
                                    nodeStates[updateNode]=UpdateMethods.updateSingleNodeBoolean(nw, pastState, updateNode);
                                    if(updateNode==input&&t>=0){nodeStates[input]=inputState;}
                                    if(updateNode==input2&&t>=0){nodeStates[input2]=inputState2;}

                                }
                            }

                            if(nodeStates[(int) indexDictionary.get("Apoptosis")]==1){Apofraction1=Apofraction1+1;}
                            if(nodeStates[(int) indexDictionary.get("Apoptosis_2")]==1){Apofraction2=Apofraction2+1;}
                            if(nodeStates[(int) indexDictionary.get("Apoptosis_3")]==1){Apofraction3=Apofraction3+1;}
                            if(nodeStates[(int) indexDictionary.get("Proliferation")]==1){Prolfraction1=Prolfraction1+1;}
                            if(nodeStates[(int) indexDictionary.get("Proliferation_2")]==1){Prolfraction2=Prolfraction2+1;}
                            if(nodeStates[(int) indexDictionary.get("Proliferation_3")]==1){Prolfraction3=Prolfraction3+1;}
                            if(nodeStates[(int) indexDictionary.get("Proliferation_4")]==1){Prolfraction4=Prolfraction4+1;}
                            if(nodeStates[Apofraction3Ind]==1){Apofraction=Apofraction+1;}else if(nodeStates[Apofraction2Ind]==1){Apofraction=Apofraction+0.5;}else if(nodeStates[Apofraction1Ind]==1){Apofraction=Apofraction+0.25;}else{}
                            if(nodeStates[Prolfraction4Ind]==1){Prolfraction=Prolfraction+1;}else if(nodeStates[Prolfraction3Ind]==1){Prolfraction=Prolfraction+0.5;}else if(nodeStates[Prolfraction2Ind]==1){Prolfraction=Prolfraction+0.25;}else if(nodeStates[Prolfraction1Ind]==1){Prolfraction=Prolfraction+0.125;}else{}

                        }

                        Apofraction1=Apofraction1/IC;
                        Apofraction2=Apofraction2/IC;
                        Apofraction3=Apofraction3/IC;
                        Prolfraction1=Prolfraction1/IC;
                        Prolfraction2=Prolfraction2/IC;
                        Prolfraction3=Prolfraction3/IC;
                        Prolfraction4=Prolfraction4/IC;
                        Apofraction=Apofraction/IC;
                        Prolfraction=Prolfraction/IC;          
                        
                        line=nw.getNames()[input]+"="+inputState+"-"+nw.getNames()[input2]+"="+inputState2+"\t"+Apofraction1+"\t"+Apofraction2+"\t"+Apofraction3+"\t"+Apofraction+"\t"+Prolfraction1+"\t"+Prolfraction2+"\t"+Prolfraction3+"\t"+Prolfraction4+"\t"+Prolfraction;
                        System.out.println(line);
                        fw.writeLine(line);
                    }
                }
            }
            
        }
        fw.close();
        
    }
            
    
}


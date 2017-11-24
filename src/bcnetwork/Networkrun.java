/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package bcnetwork;

import javax.script.ScriptException;

/**
 *
 * @author JGTZ
 * @date November 2017
 */
public class Networkrun {

    //This program simulates the ER+ breast cancer network model, starting from the
    //cancer steady state initial condition. The different java classes perform different
    //simulations:
    //NetworkSimulations - simulates the model under specified single, double, triple and no perturbations 
    //NetworkSinglePerturbationSimulations - simulates the model under all possible single node perturbations
    //NetworkDoublePerturbationSimulations - simulated the model under all possible double node perturbations
    
    /**
     * @param args the name of the TXT file where the model is. For the breast cancer model it is "BreastCancerModel.txt".
     */
    
    public static void main(String[] args) throws ScriptException {

        String filename=args[0]; //This file contains the Boolean rules of the model     
        //String filename="BreastCancerModel.txt";
        
        //Simulates the model under no perturbations and outputs the timecourse
        //of the average activity of each node
        String[] test=new String[1];
        test[0]=filename;
        NetworkSimulations.BaselineTimecourse(test);
        
        //Simulates the model in the presence of Alpelisib and outputs the timecourse
        //of the average activity of each node
        test=new String[3];
        test[0]=filename;
        test[1]="Alpelisib";test[2]="1";
        NetworkSimulations.SinglePerturbationTimecourse(test);
        
        //Simulates the model in the presence of Alpelisib and PIM, and outputs
        //the timecourse of the average activity of each node        
        test=new String[5];
        test[0]=filename;
        test[1]="Alpelisib";test[2]="1";
        test[3]="PIM";test[4]="1";
        NetworkSimulations.DoublePerturbationTimecourse(test);
        
        //Simulates the model in the presence of Alpelisib and MCL1 inhibition,
        //and outputs the timecourse of the average activity of each node        
        test=new String[5];
        test[0]=filename;
        test[1]="Alpelisib";test[2]="1";
        test[3]="MCL1";test[4]="0";
        NetworkSimulations.DoublePerturbationTimecourse(test);

        //Simulates the model in the presence of Alpelisib and Everolimus,
        //and outputs the timecourse of the average activity of each node        
        test=new String[5];
        test[0]=filename;
        test[1]="Alpelisib";test[2]="1";
        test[3]="Everolimus";test[4]="1";
        NetworkSimulations.DoublePerturbationTimecourse(test);
                
        test=new String[7];
        test[0]=filename;
        test[1]="Alpelisib";test[2]="1";
        test[3]="Palbociclib";test[4]="1";
        test[5]="Fulvestrant";test[6]="1";        
        NetworkSimulations.TriplePerturbationTimecourse(test);
        
        //Simulates the model in the presence of Alpelisib and each of 
        //every possible single node perturbation
        test=new String[1];
        test[0]=filename;                
        NetworkSinglePerturbationSimulations.main(test);
        
        //Simulates the model in the presence of Alpelisib and each of
        //every possible double node perturbation                
        NetworkDoublePerturbationSimulations.main(test);
        
    }
    

}

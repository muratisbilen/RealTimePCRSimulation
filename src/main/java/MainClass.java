import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

public class MainClass {
    public static void main(String[] args) throws Exception{
        //String refFile = "D:/Diagnovital/GISAID Data 13.01.2021 Linear.fasta";
        //String refFile = "D:/Diagnovital/cog_uk_consortium_all_sequences.fasta";
        String refFile = "GISAID_B.1.1.7_all_to_26.02.2021 Linear.fasta";

        String[] seqFiles = new String[]{"D:/Diagnovital_2/C3037T_Primer_Probe_Sequences.txt","D:/Diagnovital_2/A24403G_Primer_Probe_Sequences.txt"};

        for(int m=0;m<seqFiles.length;m++) {
            String seqFile = seqFiles[m];

            String outFile = seqFile + " - " + refFile + " Results.txt";

            BufferedReader br0 = new BufferedReader(new FileReader("D:/Diagnovital_2/"+refFile));
            BufferedReader br1 = new BufferedReader(new FileReader(seqFile));
            BufferedWriter wr = new BufferedWriter(new FileWriter(outFile));

            ArrayList<String> forw = new ArrayList<>();
            ArrayList<String> reve = new ArrayList<>();
            ArrayList<String> prob = new ArrayList<>();

            ArrayList<String> forn = new ArrayList<>();
            ArrayList<String> revn = new ArrayList<>();
            ArrayList<String> pron = new ArrayList<>();

            String s = br1.readLine();
            while ((s = br1.readLine()) != null) {
                String[] s2 = s.split("\t");
                String type = s2[1];
                if (type.equalsIgnoreCase("F")) {
                    forw.add(s2[2]);
                    forn.add(s2[0]);
                } else if (type.equalsIgnoreCase("R")) {
                    reve.add(s2[2]);
                    revn.add(s2[0]);
                } else if (type.equalsIgnoreCase("P")) {
                    prob.add(s2[2]);
                    pron.add(s2[0]);
                }
            }
            br1.close();

            String ref = "";
            String refn = "";
            String f = "";
            String r = "";
            String p = "";
            wr.write("Sample\tisDetected\tCombination\tForward Primer Pos\tReverse Primer Pos\tProbe Pos\tAmplicon Length" + "\t" + "Reason" + "\n");
            int count = 0;
            while ((refn = br0.readLine()) != null) {
                System.out.print((count++) + "\r");
                ref = br0.readLine().replace("-", "");
                for (int i = 0; i < forw.size(); i++) {
                    f = forw.get(i);
                    for (int j = 0; j < reve.size(); j++) {
                        r = reve.get(j);
                        for (int k = 0; k < prob.size(); k++) {
                            p = prob.get(k);
                            SingleSet ss = new SingleSet(forn.get(i), revn.get(j), pron.get(k), f, r, p);
                            SingleRun sr = new SingleRun(refn, ref, ss);
                            String res = refn + "\t" + sr.getDetected() + "\t" + ss.getForname() + "-" + ss.getRevname() + "-" + ss.getProname() + "\t" + sr.getPosFor() + "\t" + sr.getPosRev() + "\t" + sr.getPosPro() + "\t" + sr.getAmpliconLength() + "\t" + sr.getFailReason();
                            wr.write(res + "\n");
                        }
                    }
                }
            }
            br0.close();
            wr.close();

            BufferedReader br = new BufferedReader(new FileReader(outFile));

            s = br.readLine();
            HashMap<String, Integer> hm = new HashMap<>();
            HashMap<String, Integer> hmtotal = new HashMap<>();
            int count2 = 0;

            while ((s = br.readLine()) != null) {
                ArrayList<String> keys = new ArrayList<>(hm.keySet());
                String[] s2 = s.split("\t");
                int isd = Integer.parseInt(s2[1]);
                String comb = s2[2];
                int pos = keys.indexOf(comb);

                if (pos < 0) {
                    if (isd == 1) {
                        hm.put(comb, 1);
                    } else {
                        hm.put(comb, 0);
                    }
                    hmtotal.put(comb, 1);
                } else {
                    if (isd == 1) {
                        hm.replace(comb, hm.get(comb) + 1);
                    }
                    hmtotal.replace(comb, hmtotal.get(comb) + 1);
                }
                System.out.print((count2++) + "\r");
            }

            BufferedWriter wr2 = new BufferedWriter(new FileWriter(outFile + "_detection_rates.txt"));
            wr2.write("Combination\tDetection Rate" + "\n");
            ArrayList<String> allkeys = new ArrayList<>(hm.keySet());
            for (String ss : allkeys) {
                int d = hm.get(ss);
                int all = hmtotal.get(ss);
                wr2.write(ss + "\t" + (1.0 * d / all) + "\n");
            }
            wr2.close();
            br.close();
        }
    }

}

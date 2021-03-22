import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

public class SingleRun {
    private final int tempThreshold = 55;
    private String sequenceName = "";
    private String sequence = "";
    private SingleSet set = new SingleSet();
    private int detected = 0;
    private int ampliconLength = 0;
    private int posFor=-1, posRev=-1, posPro=-1;
    private String failReason = "";

    public SingleRun(){

    }

    public SingleRun(String sequenceName, String sequence, SingleSet set) {
        this.sequenceName = sequenceName;
        this.sequence = sequence;
        this.set = set;
        int[] res = detected(sequence,set);
        this.detected = res[0];
        this.ampliconLength = res[1];
        this.posFor = res[2];
        this.posRev = res[3];
        this.posPro = res[4];
        if(res[5]==0){
            this.failReason = "";
        }else if(res[5]==1){
            this.failReason = "Too Long Amplicon Sequence";
        }else if(res[5]==2){
            this.failReason = "No Proper Sequence Alignment";
        }else if(res[5]==3){
            this.failReason = "Sequence Overlap";
        }else{
            this.failReason = "Unknown";
        }
    }

    public int getPosFor() {
        return posFor;
    }

    public void setPosFor(int posFor) {
        this.posFor = posFor;
    }

    public int getPosRev() {
        return posRev;
    }

    public void setPosRev(int posRev) {
        this.posRev = posRev;
    }

    public int getPosPro() {
        return posPro;
    }

    public void setPosPro(int posPro) {
        this.posPro = posPro;
    }

    public String getSequenceName() {
        return sequenceName;
    }

    public void setSequenceName(String sequenceName) {
        this.sequenceName = sequenceName;
    }

    public String getSequence() {
        return sequence;
    }

    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    public SingleSet getSet() {
        return set;
    }

    public void setSet(SingleSet set) {
        this.set = set;
    }

    public int getDetected() {
        return detected;
    }

    public int getAmpliconLength() {
        return ampliconLength;
    }

    public void setAmpliconLength(int ampliconLength) {
        this.ampliconLength = ampliconLength;
    }

    public String getFailReason() {
        return failReason;
    }

    public void setFailReason(String failReason) {
        this.failReason = failReason;
    }

    public int[] detected(String sequence, SingleSet set){
        int[] res = new int[]{0,0,-1,-1,-1,0};
        String f = set.getForward().toUpperCase();
        String r = set.getReverse().toUpperCase();
        String p = set.getProbe().toUpperCase();

        int posf = findBestPosition(sequence.toUpperCase(),f,0);
        int posr = findBestPosition(sequence.toUpperCase(),r,1);
        int posp = findBestPosition(sequence.toUpperCase(),p,2);

        res[2] = posf;
        res[3] = posr;
        res[4] = posp;

        int al = posr-posf+r.length();

        res[1] = al;
        if(al>400){
            res[5] = 1; //Long Amplicon Sequence
            return(res);
        }else if(posf+f.length()-1<posp && posr>posp+p.length()-1){
            res[1] = al;
            res[0] = 1;
            return(res);
        }else if(posf<0 || posr<0 || posp<0){
            res[5] = 2; // No Sequence Alignment
            return(res);
        }else{
            res[5] = 3; // Sequence Overlap
            return(res);
        }
    }

    public int findBestPosition(String sequence, String s,int type){
        int pos = sequence.indexOf(s);
        if(pos<0){
            int[] poss = findBestPositionWithMM(sequence,s,type);
            pos = poss[0];
            int mm = poss[1];
            int gc = poss[2];
            int mmcheck = poss[3];
            double tm = getTm(sequence.substring(pos,pos+s.length()),s,(int)(100.0*mm/s.length()),(int)(100.0*gc/s.length()));
            if(mmcheck==1 && tm>tempThreshold){
                return(pos);
            }else{
                return(-1);
            }
        }
        return(pos);
    }

    public double getTm(String seq, String primer, int mm,int percGC){
        double tm = 0;
        if(primer.length() != seq.length()){
            System.out.println("Lengths of primer and seq are not the same! - getTemp()");
            return(0);
        }else {
            int psize = primer.length();

            tm = 81.5 + 0.41*percGC - 1.0*675/psize - mm;
            return (tm);
        }
    }



    public int[] findBestPositionWithMM(String sequence,String s, int type){ //Type 0: Forward, Type 1: Reverse Complement
        int[] pos = new int[]{-1,-1,-1,-1};

        String[] arr10 = new String[]{"AA","GG","CC","TT"};
        String[] arr20 = new String[]{"AG","GA","TC","CT"};

        List<String> arr1 = Arrays.asList(arr10);
        List<String> arr2 = Arrays.asList(arr20);

        int bestScore = Integer.MIN_VALUE;
        int mmcount = 0;
        int gc = 0;
        int mmcheck = 1;

        for(int i=0;i<sequence.length()-s.length();i++){
            String seq = sequence.substring(i,i+s.length());
            int score = 0;
            mmcount = 0;
            gc = 0;
            mmcheck = 1;
            for(int j=0;j<seq.length();j++){
                String x = s.substring(j,j+1);
                String y = seq.substring(j,j+1);

                int ind1 = arr1.indexOf(x+y);
                int ind2 = arr2.indexOf(x+y);

                if(ind1>-1){
                    score += 1;
                }else if(ind2>-1){
                    score += -0.5;
                }else{
                    score += -1;
                }
                if(!x.equalsIgnoreCase(y)){
                    mmcount++;
                    if(type==0 && j>seq.length()/2){
                        mmcheck=0;
                    }else if((type==1 || type==2) && j<seq.length()/2){
                        mmcheck=0;
                    }
                }
                if(x.equalsIgnoreCase("G") || x.equalsIgnoreCase("C")){
                    gc++;
                }
            }
            if(score>bestScore){
                bestScore = score;
                pos[0] = i;
                pos[1] = mmcount;
                pos[2] = gc;
                pos[3] = mmcheck;
            }
        }

        return(pos);
    }

}

import java.util.ArrayList;

public class SingleSet {
    private String forname = "";
    private String revname = "";
    private String proname = "";
    private String forward = "";
    private String reverse = "";
    private String probe = "";

    public SingleSet(){

    }
    public SingleSet(String forward, String reverse, String probe) {
        this.forward = forward;
        this.reverse = revcomp(reverse);
        this.probe = probe;
    }

    public SingleSet(String forname, String revname, String proname, String forward, String reverse, String probe) {
        this.forname = forname;
        this.revname = revname;
        this.proname = proname;
        this.forward = forward;
        this.reverse = revcomp(reverse);
        this.probe = probe;
    }

    public String getForward() {
        return forward;
    }

    public void setForward(String forward) {
        this.forward = forward;
    }

    public String getReverse() {
        return reverse;
    }

    public void setReverse(String reverse) {
        this.reverse = revcomp(reverse);
    }

    public String getProbe() {
        return probe;
    }

    public void setProbe(String probe) {
        this.probe = probe;
    }

    public String getForname() {
        return forname;
    }

    public void setForname(String forname) {
        this.forname = forname;
    }

    public String getRevname() {
        return revname;
    }

    public void setRevname(String revname) {
        this.revname = revname;
    }

    public String getProname() {
        return proname;
    }

    public void setProname(String proname) {
        this.proname = proname;
    }

    public static String revcomp(String x){
        String x2 = "";
        int len = x.length();
        String[] lets = new String[]{"A","C","G","T","U","W","S","M","K","R","Y","B","D","H","V","N","Z"};
        String[] lets2 = new String[]{"T", "G","C","A","A","W","S","K","M","Y","R","V","H","D","B","N","Z"};
        ArrayList let1 = new ArrayList<>();
        ArrayList let2 = new ArrayList<>();

        for(int i=0;i<lets.length;i++){
            let1.add(lets[i]);
            let2.add(lets2[i]);
        }

        for(int i=len;i>0;i--){
            String a = x.substring(i-1,i);
            int pos = let1.indexOf(a);
            if(pos>-1){
                x2 += let2.get(pos);
            }else{
                x2 += a;
            }
        }
        return(x2);
    }

}

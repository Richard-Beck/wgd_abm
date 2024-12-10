package metastasis;

import HAL.Rand;

import java.util.ArrayList;
import java.util.Iterator;

// rationale for this class is that we have implemented cells as a 2d grid agent but this
// ceases to make sense after intravasation.
public class CTCs {
    int wgd;
    double pDie=0.0;

    ArrayList<Boolean> ctcs;

    public CTCs(){
        ctcs = new ArrayList<Boolean>();
    }

    public void stepCTCs(Rand rn){
        Iterator<Boolean> iterator = ctcs.iterator();
        while (iterator.hasNext()) {
            Boolean value = iterator.next();
            if (rn.Double()<pDie) {
                iterator.remove();
            }
        }
    }


}

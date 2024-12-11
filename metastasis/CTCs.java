package metastasis;

import HAL.Rand;

import java.util.ArrayList;
import java.util.Iterator;

// rationale for this class is that we have implemented cells as a 2d grid agent but this
// ceases to make sense after intravasation.
public class CTCs {
    double pDie=0.01;
    double pEx = 0.001;

    ArrayList<CellEx> ctcList;

    public CTCs(){
        ctcList = new ArrayList<>();
    }

    public ArrayList<CellEx> stepCTCs(Rand rn){
        ArrayList<CellEx> extravasators = new ArrayList<>();
        Iterator<CellEx> iterator = ctcList.iterator();
        while (iterator.hasNext()) {
            // does this line do anything?
            CellEx value = iterator.next();
            if (rn.Double()<pDie) {
                iterator.remove();
            }else if(rn.Double()<pEx){
                extravasators.add(value);
                iterator.remove();
            }

        }
        return extravasators;
    }
}

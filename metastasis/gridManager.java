package metastasis;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.UIGrid;
import HAL.Util;

public class gridManager {
    PDEGrid2D MMP;
    PDEGrid2D ECM;

    int size;

    double DegMMP = 0.1;
    double MMPDiffCoef=0.1;
    double DegECMbyMMP=0.1;

    public gridManager(int x, int y) {
        MMP = new PDEGrid2D(x, y);
        ECM = new PDEGrid2D(x, y);
        size = x*y;
        for(int i = 0; i<size; i++){
            ECM.Set(i,1.);
        }
        ECM.Update();
    }

    void diffuse(){
        MMP.DiffusionADI(MMPDiffCoef);
        MMP.Update();
        for(int i = 0; i< size;i++){
            MMP.Add(i,-MMP.Get(i)*DegMMP);
        }
        MMP.Update();
        ECM.Update();
    }

    void degradeECM(){
        // this models only the degredation of ECM via MMP
        // cells also degrade ECM, handled elsewhere.
        for(int i = 0; i<size; i++){

            ECM.Add(i,-ECM.Get(i)*Math.min(1,MMP.Get(i)*DegECMbyMMP));
        }
        ECM.Update();
    }


}

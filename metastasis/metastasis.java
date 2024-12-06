package metastasis;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.AgentSQ2D;
import HAL.Gui.GridWindow;
import HAL.Gui.UIGrid;
import HAL.Tools.FileIO;
import HAL.Rand;
import HAL.Tools.PhylogenyTracker.Genome;
import HAL.Util;

import java.util.*;

//cells grow and mutate
class CellEx extends AgentSQ2D<metastasis>{
    int nMutations;
    int wgd = 0;
    double MMPGen=.01;
    double ECMDeg=0.01;
    int popMax=4; // max number of agents that can populate a site

    double ECMGradX,ECMGradY;
    double HaptoCoef=0.001;

    double pMove = 0.001;
    double pwgd = 0.001;
    void Mutate(){
        if(nMutations< G.MAX_MUTATIONS && G.rn.Double()< G.MUT_PROB){
            nMutations++;
        }
        if(G.rn.Double()<pwgd){
            wgd=1-wgd;
            if(wgd==1){
                //pMove=pMove*2;
                HaptoCoef=HaptoCoef*2;
            }else{
                //pMove=pMove/2;
                HaptoCoef=HaptoCoef/2;
            }
        }
    }

    void Draw(int pop){
        //double shade = (double) pop /popMax;
        //G.vis.SetPix(Xsq(), Ysq(), Util.RGB(shade, shade,shade));//sets a single pixel
    }

    void Divide(){
        int pop = G.PopAt(Isq());//finds von neumann neighborhood indices around cell.
        if(pop<popMax){
            CellEx daughter= G.NewAgentSQ(Isq());//generate a daughter, the other is technically the original cell
            daughter.nMutations=nMutations;//start both daughters with same number of mutations
            daughter.Draw(pop);
            Mutate();//during division, there is a possibility of mutation of one or both daughters
            daughter.Mutate();
        }
    }


    void interactWithFields(gridManager M){
        M.ECM.Add(Isq(),-ECMDeg*M.ECM.Get(Isq()));
        M.MMP.Add(Isq(), MMPGen);
        ECMGradX = M.ECM.Get(Math.min(Xsq()+1,G.xDim-1),Ysq())-M.ECM.Get(Math.max(Xsq()-1,0),Ysq());
        ECMGradY = M.ECM.Get(Xsq(),Math.min(Ysq()+1,G.yDim-1))-M.ECM.Get(Xsq(),Math.max(Ysq()-1,0));
    }
    int selectIndex(double[] probs) {
        double r = G.rn.Double() * Arrays.stream(probs).sum(), cum = 0;
        for (int i = 0; i < probs.length; i++) if ((cum += probs[i]) >= r) return i;
        return -1;
    }
    void tryMove(){
        // try move due to diffusion
        double[] pq={Math.max(0,ECMGradY*HaptoCoef)+pMove/4,
                Math.max(0,ECMGradX*HaptoCoef)+pMove/4,
                Math.max(0,-ECMGradX*HaptoCoef)+pMove/4,
                Math.max(0,-ECMGradY*HaptoCoef)+pMove/4,
                1-pMove};
        double pTotal = 0.;
        for(double d :pq) pTotal+=d;
        for(int i = 0; i<pq.length;i++) pq[i]/=pTotal;

        int newIndex = selectIndex(pq);

        if(newIndex==0){
            MoveSafeSQ(G.ItoX(Isq()),G.ItoY(Isq())+1);
        }
        if(newIndex==1){
            MoveSafeSQ(G.ItoX(Isq())+1,G.ItoY(Isq()));
        }
        if(newIndex==2){
            MoveSafeSQ(G.ItoX(Isq()),G.ItoY(Isq())-1);
        }
        if(newIndex==3){
            MoveSafeSQ(G.ItoX(Isq())-1,G.ItoY(Isq()));
        }
    }
}

class vessels {
    List<int[]> coords;
    double vesselValue = 10.;

    public vessels(int N, int xmax, int ymax, Rand rn) {
        coords = new ArrayList<>();
        for (int i = 0; i < N; i++) {

            int x = rn.Int(xmax);
            int y = rn.Int(ymax);
            coords.add(new int[]{x, y});
        }
    }

    //void interactWithFields(gridManager M){
      //  for(int[] coord : coords){
        //    M.O2.Set(coord[0],coord[1],vesselValue);
        //}
    //}
}

public class metastasis extends AgentGrid2D<CellEx> {
    final static int BLACK= Util.RGB(0,0,0);
    double DIV_PROB =0.2;
    double MUT_PROB =0.003;
    double DIE_PROB =0.05;
    double MUT_ADVANTAGE =1.08;
    int MAX_MUTATIONS =19;
    int[]mutCounts=new int[MAX_MUTATIONS+1];//+1 to count for un-mutated type
    int[]hood=Util.GenHood2D(new int[]{1,0,-1,0,0,1,0,-1}); //equivalent to int[]hood=Util.VonNeumannHood(false);
    Rand rn=new Rand(1);
    UIGrid vis;
    FileIO outputFile=null;

    gridManager M;

    public metastasis(int x, int y, UIGrid vis) {
        super(x, y, CellEx.class);
        this.vis=vis;
    }
    public metastasis(int x, int y, UIGrid vis, String outputFileName) {
        super(x, y, CellEx.class);
        this.vis=vis;
        outputFile=new FileIO(outputFileName,"w");
    }
    public void InitTumor(double radius){
        //places tumor cells in a circle
        int[]circleHood= Util.CircleHood(true,radius);//generate circle neighborhood [x1,y1,x2,y2,...]
        int len=MapHood(circleHood,xDim/2,yDim/2);
        for (int i = 0; i < len; i++) {
            CellEx c=NewAgentSQ(circleHood[i]);
            c.nMutations=0;
            c.Draw(1);
        }
    }

    public void drawAllCells(){
        Map<List<Integer>, double[]> map = new HashMap<>();
        for (CellEx c : this) {
            List<Integer> key = Arrays.asList(c.Xsq(), c.Ysq());
            double[] value = new double[]{(double) c.wgd / c.popMax, (double)(1 - c.wgd) / c.popMax, (double) c.wgd / c.popMax};

            map.merge(key, value, (oldValue, newValue) -> {
                for (int i = 0; i < oldValue.length; i++) {
                    oldValue[i] += newValue[i];
                }
                return oldValue;
            });
        }
        for (Map.Entry<List<Integer>, double[]> entry : map.entrySet()) {
            List<Integer> key = entry.getKey();
            double[] value = entry.getValue();
            vis.SetPix(key.get(0), key.get(1), Util.RGB(value[0], value[1], value[2]));
        }


    }
    public void StepCells(int tick){
        Arrays.fill(mutCounts,0);//clear the mutation counts
        for (CellEx c : this) {//iterate over all cells in the grid
            mutCounts[c.nMutations]++;//count up all cell types for this timestep
            c.interactWithFields(M);
            c.tryMove();
            if(rn.Double()< DIE_PROB){
                int newPop = PopAt(c.Isq())-1;
                c.Draw(newPop);
                c.Dispose();//removes cell from sptial grid and iteration
            }
            else if(rn.Double()< DIV_PROB*Math.pow(MUT_ADVANTAGE,c.nMutations)){//application of mutational advantage
                c.Divide();
            }
        }
        if(outputFile!=null){
            outputFile.Write(Util.ArrToString(mutCounts,",")+"\n");//write populations every timestep
        }
        ShuffleAgents(rn);//shuffles order of for loop iteration
//        IncTick();//increments timestep, including newly generated cells in the next round of iteration
    }

    public void drawMMP(int xshift,int yshift){
        for(int i =0; i<(xDim*yDim); i++){
            vis.SetPix(ItoX(i)+xshift,ItoY(i)+yshift,Util.RGB(M.MMP.Get(i),0,0));
        }
    }
    public void drawECM(int xshift,int yshift){
        for(int i =0; i<(xDim*yDim); i++){
            vis.SetPix(ItoX(i)+xshift,ItoY(i)+yshift,Util.RGB(0,0,M.ECM.Get(i)));
        }
    }
    public static void main(String[]args){
        ArrayList<Double[]>out=new ArrayList<>();
        //int x=500,y=500,scaleFactor=2;
        int x=200,y=200,scaleFactor=2;
        GridWindow vis=new GridWindow(3*x,y,scaleFactor);//used for visualization
        metastasis grid=new metastasis(x,y,vis);
        grid.InitTumor(5);
        grid.M = new gridManager(x,y);
        vessels V = new vessels(100,x,y, grid.rn);

        for (int tick = 0; tick < 20000; tick++) {
            vis.TickPause(0);//set to nonzero value to cap tick rate.
            grid.drawMMP(x,0);
            grid.drawECM(2*x,0);
            grid.drawAllCells();
            grid.StepCells(tick);
            grid.M.degradeECM();
            grid.M.diffuse();

        }
    }
}
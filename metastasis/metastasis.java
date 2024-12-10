package metastasis;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2D;
import HAL.Gui.OpenGL2DWindow;
import HAL.Tools.FileIO;
import HAL.Rand;
import HAL.Util;

import java.util.*;

//cells grow and mutate
class CellEx extends AgentSQ2D<metastasis>{
    int nMutations;
    double ECMGradX,ECMGradY;
    double pwgd = 0.01;

    double pIn = 0.001;

    public enum CellType {

        tumor, wgd, vessel;
        double MMPGen,ECMDeg,HaptoCoef,pMove,space;
        int color;
        // Constructor for initializing the parameters
        public void setParams(double MMPGen, double ECMDeg, double HaptoCoef,double pMove,int color, double space) {
            this.color=color;
            this.MMPGen=MMPGen;
            this.ECMDeg=ECMDeg;
            this.HaptoCoef=HaptoCoef;
            this.pMove=pMove;
            this.space=space;
        }

    }

    CellType type;
    public static void initialize_types(){
        CellType.valueOf("wgd").setParams(0.01,0.01,0.002,0.01,Util.RGB(0,1,0),0.5);
        CellType.valueOf("tumor").setParams(0.01,0.01,0.002,0.002,Util.RGB(0,0,1),0.25);
        CellType.valueOf("vessel").setParams(0.0,0.0,0.0,0.0,Util.RGB(1,0,0),1);
    }

    public void Init(CellType type) {
        this.type = type;
        // Any parameters that belong to the base clase could be set below:
    }

    void Divide(){
        ArrayList<CellEx> agents = new ArrayList<>();
        G.GetAgents(agents,Isq());
        double space = 0;
        for(CellEx c:agents) space=space+c.type.space;
        if((1-space)>=type.space){
            if(type==CellType.wgd | G.rn.Double()>pwgd){
                CellEx daughter= G.NewAgentSQ(Isq());//generate a daughter, the other is technically the original cell
                daughter.Init(type);
            } else{
                Init(CellType.wgd);
            }
        }
    }

    boolean tryIntravasate(){
        // if a vessel can share a site with a tumor cell we'd need to use getAgents()
        if(G.rn.Double()< pIn){
            int nopts = MapHood(G.hood);
            int trial = G.hood[G.rn.Int(nopts)];
            if(G.PopAt(trial)==0) return false;
            if(G.GetAgent(trial).type==CellType.vessel){
                return true;
            }
        }
        return false;
    }

    void interactWithFields(gridManager M){
        M.ECM.Add(Isq(),-type.ECMDeg*M.ECM.Get(Isq()));
        M.MMP.Add(Isq(), type.MMPGen);
        // IDK if its best for the cells to know these things or not. Think on it again when code is more done
        ECMGradX = M.ECM.Get(Math.min(Xsq()+1,G.xDim-1),Ysq())-M.ECM.Get(Math.max(Xsq()-1,0),Ysq());
        ECMGradY = M.ECM.Get(Xsq(),Math.min(Ysq()+1,G.yDim-1))-M.ECM.Get(Xsq(),Math.max(Ysq()-1,0));
    }
    int selectIndex(double[] probs) {
        double r = G.rn.Double() * Arrays.stream(probs).sum(), cum = 0;
        for (int i = 0; i < probs.length; i++) if ((cum += probs[i]) >= r) return i;
        return -1;
    }
    void tryMove(){
        // try move due to diffusion & haptotaxis
        // must be a better way to write this
        double[] pq={Math.max(0,ECMGradY*type.HaptoCoef)+type.pMove/4,
                Math.max(0,ECMGradX*type.HaptoCoef)+type.pMove/4,
                Math.max(0,-ECMGradX*type.HaptoCoef)+type.pMove/4,
                Math.max(0,-ECMGradY*type.HaptoCoef)+type.pMove/4,
                1-type.pMove};
        double pTotal = 0.;
        for(double d :pq) pTotal+=d;
        for(int i = 0; i<pq.length;i++) pq[i]/=pTotal;

        int newIndex = selectIndex(pq);

        int newX = Xsq();
        int newY = Ysq();
        switch (newIndex) {
            case 0:
                newY = newY+1; break;
            case 1:
                newX = newX+1; break;
            case 2:
                newY = newY-1; break;
            case 3:
                newX = newX-1; break;
            default:
                break;
        }

        if(G.In(newX,newY)){
            ArrayList<CellEx> agents = new ArrayList<>();
            G.GetAgents(agents,newX,newY);
            double space = 0;
            for(CellEx c:agents) space=space+c.type.space;
            if(1-space>=type.space) MoveSQ(newX,newY);
        }
    }
}

public class metastasis extends AgentGrid2D<CellEx> {
    final static int BLACK= Util.RGB(0,0,0);
    double DIV_PROB =0.2;
    double DIE_PROB =0.02;
    double MUT_ADVANTAGE =1.08;
    int MAX_MUTATIONS =19;
    int[]mutCounts=new int[MAX_MUTATIONS+1];//+1 to count for un-mutated type
    int[]hood=Util.GenHood2D(new int[]{1,0,-1,0,0,1,0,-1}); //equivalent to int[]hood=Util.VonNeumannHood(false);
    Rand rn=new Rand(1);
    OpenGL2DWindow vis;
    FileIO outputFile=null;

    gridManager M;
    CTCs ctcs=new CTCs();
    public metastasis(int x, int y, OpenGL2DWindow vis) {
        super(x, y, CellEx.class);
        this.vis=vis;
    }
    public metastasis(int x, int y, OpenGL2DWindow vis, String outputFileName) {
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
            c.Init(CellEx.CellType.tumor);
        }
        // also initialize the vessels here?
        int N = 100;
        for (int i = 0; i < N; i++) {
            int x = rn.Int(xDim);
            int y = rn.Int(yDim);
            CellEx c=NewAgentSQ(x,y);
            c.Init(CellEx.CellType.vessel);
        }
    }

    public void StepCells(int tick){
        for (CellEx c : this) {//iterate over all cells in the grid
            if(c.type== CellEx.CellType.vessel) {
                // ... vessel stuff
                continue;
            }
            if(c.tryIntravasate()){
                ctcs.ctcs.add(c.type== CellEx.CellType.wgd);
                c.Dispose();
                continue;
            }
            c.interactWithFields(M);
            c.tryMove();
            if(rn.Double()< DIE_PROB){
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
    // these are utility functions for adding colors, perhaps should be elsewhere.
    public static int AddARGB(int color1, int color2) {
        int r = Math.min(255, ((color1 >> 16) & 0xff) + ((color2 >> 16) & 0xff));
        int g = Math.min(255, ((color1 >> 8) & 0xff) + ((color2 >> 8) & 0xff));
        int b = Math.min(255, (color1 & 0xff) + (color2 & 0xff));
        return 0xff000000 | (r << 16) | (g << 8) | b;
    }

    public static int MulARGB(int color, double factor) {
        int r = Math.min(255, (int)(((color >> 16) & 0xff) * factor));
        int g = Math.min(255, (int)(((color >> 8) & 0xff) * factor));
        int b = Math.min(255, (int)((color & 0xff) * factor));
        return 0xff000000 | (r << 16) | (g << 8) | b;
    }
    public void draw(){
        for(int i = 0; i< length; i++){
            ArrayList<CellEx> agents = new ArrayList<>();
            GetAgents(agents,i);
            int value = BLACK;
            for(CellEx a : agents){
                value = AddARGB(value,MulARGB(a.type.color,a.type.space));
            }
            vis.SetPix(i,value);
            vis.SetPix(i + length,Util.RGB(M.MMP.Get(i),0,0));
            vis.SetPix(i + 2*length,Util.RGB(0,0,M.ECM.Get(i)));
        }
        vis.Update();
    }

    public static void main(String[]args){
        ArrayList<Double[]>out=new ArrayList<>();
        //int x=500,y=500,scaleFactor=2;
        int x=200,y=200,scaleFactor=2;
        OpenGL2DWindow vis = new OpenGL2DWindow("Primary", 9*x, 3*y, 3*x, y);
        metastasis grid=new metastasis(x,y,vis);
        CellEx.initialize_types();
        grid.InitTumor(5);
        grid.M = new gridManager(x,y);
        for (int tick = 0; tick < 20000; tick++) {
            vis.TickPause(0);//set to nonzero value to cap tick rate.
            grid.draw();
            grid.ctcs.stepCTCs(grid.rn);
            System.out.println(grid.ctcs.ctcs.size());
            grid.StepCells(tick);
            grid.M.degradeECM();
            grid.M.diffuse();
            if(vis.IsClosed()) break;
        }
    }
}
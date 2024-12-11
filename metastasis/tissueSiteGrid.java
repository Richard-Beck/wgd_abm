package metastasis;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2D;
import HAL.Tools.FileIO;
import HAL.Rand;
import HAL.Util;

import java.util.*;
import java.util.stream.Collectors;

class CellEx extends AgentSQ2D<tissueSiteGrid>{
    int nMutations;
    double ECMGradX,ECMGradY;
    double pwgd = 0.01;

    double pIn = 0.001;

    BitSet genome;

    double pMut=0.1;

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

    public void Init(CellType type,BitSet genome) {
        this.type = type;
        this.genome = genome;
        // Any parameters that belong to the base clase could be set below:
    }

    int Divide(){
        ArrayList<CellEx> agents = new ArrayList<>();
        G.GetAgents(agents,Isq());
        double space = 0;
        int nNew=0;
        for(CellEx c:agents) space=space+c.type.space;
        if((1-space)>=type.space){
            if(type==CellType.wgd | G.rn.Double()>pwgd){
                CellEx daughter= G.NewAgentSQ(Isq());//generate a daughter, the other is technically the original cell
                daughter.Init(type,genome);
                nNew++;
            } else{
                Init(CellType.wgd,genome);
            }
        }
        return nNew;
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
        M.ECM.Add(Isq(),-type.space*type.ECMDeg*M.ECM.Get(Isq()));
        M.MMP.Add(Isq(), type.space*type.MMPGen);
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
public class tissueSiteGrid extends AgentGrid2D<CellEx> {
    final static int BLACK= Util.RGB(0,0,0);
    double DIV_PROB =0.2;
    double pExtr=0.01;
    double DIE_PROB =0.02;
    double MUT_ADVANTAGE =1.08;
    int MAX_MUTATIONS =19;

    int tumorPop;
    int[]mutCounts=new int[MAX_MUTATIONS+1];//+1 to count for un-mutated type
    int[]hood=Util.GenHood2D(new int[]{1,0,-1,0,0,1,0,-1}); //equivalent to int[]hood=Util.VonNeumannHood(false);
    Rand rn;
    Viz vis;
    FileIO outputFile=null;
    gridManager M;

    public tissueSiteGrid(int x, int y, Viz vis, Rand globalRNG) {
        super(x, y, CellEx.class);
        this.rn=globalRNG;
        this.vis=vis;
    }
    public tissueSiteGrid(int x, int y, Viz vis, String outputFileName) {
        super(x, y, CellEx.class);
        this.vis=vis;
        outputFile=new FileIO(outputFileName,"w");
    }
    public void InitTumor(double radius){
        //places tumor cells in a circle
        if(radius>0){
            int[]circleHood= Util.CircleHood(true,radius);//generate circle neighborhood [x1,y1,x2,y2,...]
            int len=MapHood(circleHood,xDim/2,yDim/2);
            for (int i = 0; i < len; i++) {
                CellEx c=NewAgentSQ(circleHood[i]);
                BitSet genome = new BitSet(100);
                c.Init(CellEx.CellType.tumor,genome);
            }
        }

        // also initialize the vessels here?
        int N = 100;
        for (int i = 0; i < N; i++) {
            int x = rn.Int(xDim);
            int y = rn.Int(yDim);
            CellEx c=NewAgentSQ(x,y);
            BitSet genome = new BitSet(100);
            c.Init(CellEx.CellType.vessel,genome);
        }
    }

    public void StepCells(int tick,CTCs C){
        tumorPop=0;
        for (CellEx c : this) {//iterate over all cells in the grid
            if(c.type== CellEx.CellType.vessel) {
               // if(rn.Double()<pExtr*C.ctcList.size()){
                 //   int i = rn.Int(C.ctcList.size());
                   // CellEx e = NewAgentPT(c.Xsq(),c.Ysq());
                    //e.Init(C.ctcList.get(i).type,C.ctcList.get(i).genome);
                    //C.ctcList.remove(i);
                //}
                // ... vessel stuff
                continue;
            }
            tumorPop++;
            if(rn.Double()<c.pMut) {
                c.genome.flip(rn.Int(c.genome.size()));
            }
            if(c.tryIntravasate()){
                C.ctcList.add(c);
                c.Dispose();
                tumorPop--;
                continue;
            }
            c.interactWithFields(M);
            c.tryMove();
            if(rn.Double()< DIE_PROB){
                c.Dispose();//removes cell from sptial grid and iteration
                tumorPop--;
            }
            else if(rn.Double()< DIV_PROB*Math.pow(MUT_ADVANTAGE,c.nMutations)){//application of mutational advantage
                tumorPop+=c.Divide();
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
    public int draw(){
        vis.makeContextCurrent();
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
        return vis.Update();
    }

    public void writeData(FileIO writer,int gridNo){
        long[] genomeArray;
        for(CellEx c:this){
            int x = c.Xsq();
            int y = c.Ysq();
            genomeArray = c.genome.toLongArray();
            String genomeString = Arrays.stream(genomeArray)
                    .mapToObj(String::valueOf)
                    .collect(Collectors.joining("|"));

            writer.Write(String.valueOf(c.Xsq())+","+
                    String.valueOf(c.Ysq())+","+
                    genomeString+","+
                    c.type+","+gridNo+"\n");
        }

    }

    public void stepTumor(int tick, CTCs C,boolean drawMe){
        if(tick%20==0 & drawMe) {
            draw();
        }
        StepCells(tick, C);
        M.degradeECM();
        M.diffuse();
    }
}
package Examples._3OffLatticeExample;

import HAL.GridsAndAgents.PDEGrid2D;
import HAL.GridsAndAgents.SphericalAgent2D;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import HAL.Gui.OpenGL2DWindow;
import HAL.Tools.FileIO;
import HAL.Tools.Internal.Gaussian;
import HAL.Rand;

import java.io.IOException;
import java.util.ArrayList;

import static HAL.Util.*;


class CellOL extends SphericalAgent2D<CellOL,ExampleOffLattice>{

    double forceSum;//used with contact inhibition calculation

    public enum CellType {
        //normal(0.25, 0.8, 0.02,RGB256(0,250,250),0.2),
       // tumor(0.25, 0.05, 0.02,RGB256(200,200,100),0.5),
       // wgd(0.35, 0.05, 0.02,RGB256(200,0,200),0.5),
       // vessel(0.5, 2.5, 20.0,RGB256(255,2,3),0.01),
       // immune(0.15, 0.01, 0.005,RGB256(0,255,0),0.9);
        normal, tumor, wgd, vessel, immune;
        // Final fields for the parameters
        double radius, inhib_weight, div_bias, friction;
        int color;
        // Constructor for initializing the parameters
        public void setParams(double radius, double inhib_weight, double div_bias, int color,double friction) {
            this.radius = radius;
            this.inhib_weight = inhib_weight;
            this.div_bias = div_bias;
            this.color=color;
            this.friction=friction;
        }

    }

    // Cell type (th)
    CellType type;
    double div_bias;

    // Constructor to set the type
    public void Init(CellType type) {
        this.type = type;
        // the following is necessary because radius is part of the base class and is used in
        // some of the base methods:
        this.radius = type.radius;
        this.div_bias=type.div_bias;
    }

    double ForceCalc(double overlap,CellOL other){
        if(overlap<0) {
            return 0;//if cells aren't actually overlapping, then there is no force response
        }
        return G.FORCE_SCALER*overlap;//this constant scaling of the overlap is called Hooke's law!
    }
    public void CalcMove(){
        //sets x and y velocity components of cell
        //G.neighborList.clear();
        //G.GetAgentsRad(G.neighborList,G.neighborInfo,Xpt(),Ypt(),G.RADIUS*2);
        //NOTE - in the below it is assumed VESSELS have the largest possible radius
        forceSum=SumForces(radius+ CellType.vessel.radius,this::ForceCalc);
        CapVelocity(1.0);
        //forceSum=SumForces(G.neighborList,G.neighborInfo,this::ForceCalc);
    }

    public void CalcStim(){
        G.neighborList.clear();
        G.GetAgentsRad(G.neighborList,G.neighborInfo,Xpt(),Ypt(),radius+ CellType.wgd.radius);
        boolean is_stimulated = false;
        for (CellOL a : G.neighborList) {
            if (a != this) {
                double xComp = DispX(a.Xpt());
                double yComp = DispY(a.Ypt());
                double dist = Norm(xComp, yComp);
                double touchDist = (radius + a.radius) - dist;
                if(touchDist>0 & (a.type==CellType.wgd | a.type==CellType.tumor)){
                    is_stimulated=true;
                    if(G.rn.Double()<0.2) a.Dispose(); // KILLING!
                }

            }
        }
        if(is_stimulated) div_bias= type.div_bias;
        if(!is_stimulated) div_bias=0;
    }
    public boolean CanDivide(double div_bias,double inhib_weight){
        return G.rn.Double()<Math.tanh(div_bias-forceSum*inhib_weight);
    }
    public void MoveDiv(double Rconc){
        //move cell and reduce x and y velocity components by friction constant

        ApplyFriction(type.friction);
        ForceMove();
        if(type==CellType.vessel) return;
        if(type==CellType.immune){
            CalcStim();
        }
        if(G.rn.Double() < 0.002) Dispose();
        //compute whether division can occur, using the constants
        if(CanDivide(div_bias, type.inhib_weight)){
            if(Rconc> G.resource_threshold){
                Divide(radius*2.0/3.0, G.divCoordStorage, G.rn).Init(type);
            }
        }
    }

    public void modify_pde_grid() {
        switch (type) {
            case vessel:
                G.R.Set(Isq(), 10.);
                break;
            default:
                G.R.Mul(Isq(), -.01);
                break;
        }
    }

    public static void initializeCellTypes(String filePath)  {
        FileIO reader = new FileIO(filePath,"r");
        String line = reader.ReadLine(); // Skip header

        while ((line = reader.ReadLine()) != null) {
            String[] values = line.split(",");
            String type = values[0];

            double radius = Double.parseDouble(values[1]);
            double inhib_weight = Double.parseDouble(values[2]);
            double div_bias= Double.parseDouble(values[3]);
            int r = Integer.parseInt(values[4]);
            int g = Integer.parseInt(values[5]);
            int b = Integer.parseInt(values[6]);
            double friction = Double.parseDouble(values[7]);

            int color = RGB256(r, g, b);

            // Initialize the corresponding enum value
            CellType.valueOf(type).setParams(radius, inhib_weight, div_bias, color, friction);

        }

        reader.Close();
    }
}

public class ExampleOffLattice extends AgentGrid2D<CellOL> {


    static final int WHITE=RGB256(248,255,252), CYTOPLASM=RGB256(200,210,210),
            normal=RGB256(0,250,250),
    tumor=RGB256(200,200,100),
    wgd=RGB256(200,0,200),
    vessel=RGB256(255,2,3);;



    double RADIUS=0.25;



    double FORCE_SCALER=0.25;//this constant was found to be rather stable, but tweak it and see what happens!
    double resource_threshold = .1;

    ArrayList<CellOL> neighborList=new ArrayList<>();
    ArrayList<double[]> neighborInfo=new ArrayList<>();
    double[]divCoordStorage=new double[2];

    Rand rn=new Rand(System.currentTimeMillis());
    Gaussian gn =new Gaussian();

    public PDEGrid2D R;


    public ExampleOffLattice(int x, int y) {
        super(x, y, CellOL.class,true,true);
        R = new PDEGrid2D(xDim, yDim);
    }
    public static void main(String[] args) throws IOException {
       // OpenGL2DWindow.MakeMacCompatible(args);
        String path2params = args[0];
        String outDir = args[1];
        CellOL.initializeCellTypes(path2params);
        int x=40,y=20;
        //to record output, call the constructor with an output filename
        ExampleOffLattice ex=new ExampleOffLattice(x,y);
        ex.wrapX=false;

        OpenGL2DWindow vis = null;
        GifMaker gif = null;
        boolean immune_inf = false;
        boolean make_gif=false;
        IO writer = new IO();
       // gif.AddFrame(vis);
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-v")) {
                vis = new OpenGL2DWindow("Off Lattice Example", 25 * x, 2 * 25 * y, x, y * 2);
            }
            if (args[i].equals("-i")) {
                immune_inf = true;
            }
            if (args[i].equals("-o")) {
                writer = new IO(outDir,path2params);
            }
        }

        ex.Setup(1000,y,x,CellOL.CellType.normal);
        int init_vessels = 40;
        ex.Setup(init_vessels,y,x,CellOL.CellType.vessel);
        int i=0;
        while(i<200 && (vis == null || !vis.IsClosed())) {//check for click on close button on window
            if (vis != null) {
                vis.TickPause(0);
            }
            ex.StepCells();
            if(immune_inf) ex.immuneIn();
            ex.DrawModel(vis);
            writer.Draw(i,vis);
            writer.RecordOut(ex,i);
            i++;
        }
        ex.Setup(50,y,0.1*x,CellOL.CellType.wgd);
        ex.Setup(50,y,0.1*x,CellOL.CellType.tumor);

        while(i<10000 && (vis == null || !vis.IsClosed())){//check for click on close button on window
            if (vis != null) {
                vis.TickPause(0);
            }
            ex.StepCells();
            if(immune_inf) ex.immuneIn();
            ex.DrawModel(vis);
            if(vis!=null & i%500==0){
                //vis.ToJPG("Examples/_3OffLatticeExample/img/"+i+".jpg");
            }
            double front = ex.calcFront();
            if(front>.6*x){
                System.out.println(front);
                ex.shiftFront(x,y,0.8);
                front = ex.calcFront();
                System.out.println(front);
                // add new vessels to replace those lost. Need to put this somewhere else more organised.
                double vessel_density= (1-0.8)*init_vessels;
                int n_add = (int) Math.round(vessel_density);

                for (int ix = 0; ix < n_add; ix++) {
                    double vx = 0.8*x+ex.rn.Double(x*(1-0.8));
                    double vy = ex.rn.Double(y);
                    ex.NewAgentPT(vx,vy).Init(CellOL.CellType.vessel);
                }
            }
            writer.Draw(i,vis);
            writer.RecordOut(ex,i);
            i++;
        }

        writer.Close();
        if (vis != null) {
            vis.Close();
        }

    }

    public void Setup(double initPop, double ymax, double xmax, CellOL.CellType type){
        for (int i = 0; i < initPop; i++) {
            double x = rn.Double(xmax);
            double y = rn.Double(ymax);
            //create a new agent, and set the type depending on a comparison with the random number generator
            NewAgentPT(x,y).Init(type);
        }
    }

    public void immuneIn(){

        for(CellOL cell:this){
            if(cell.type== CellOL.CellType.vessel){
                if(rn.Double() < 0.002){
                    cell.Divide(cell.radius*2.0/3.0, divCoordStorage, rn).Init(CellOL.CellType.immune);
                };
            }
        }
    }

    public void DrawModel(OpenGL2DWindow vis){
        if(vis==null) return;
        //vis.Clear(WHITE);
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                vis.SetPix(x, y+yDim, HeatMapRGB(R.Get(x, y)));
                vis.SetPix(x, y, WHITE);

            }
        }
        for (CellOL cell : this) {
            //draw "cytoplasm" of cell
            vis.Circle(cell.Xpt(),cell.Ypt(),cell.type.radius,CYTOPLASM);
        }
        for (CellOL cell : this) {
            //draw colored "nucleus" on top of cytoplasm
            vis.Circle(cell.Xpt(), cell.Ypt(), cell.type.radius / 3, cell.type.color);
        }

        vis.Update();
    }
    public void StepCells(){
        for (CellOL cell : this) {
            cell.CalcMove();//calculation of forces before any agents move, for simultaneous movement
        }
        for (CellOL cell : this) {
            double Rconc = this.R.Get(cell.Isq());
            cell.MoveDiv(Rconc);//movement and division
            cell.modify_pde_grid();
        }
        this.R.DiffusionADI(.1);
        this.R.Update();
    }

    public double calcFront(){
        double xmax=0;
        for (CellOL cell : this) {
            if(cell.Xpt()>xmax & cell.type==CellOL.CellType.tumor) xmax = cell.Xpt();
            if(cell.Xpt()>xmax & cell.type==CellOL.CellType.wgd) xmax = cell.Xpt();
        }
        return(xmax);
    }

    public double calcFrontDEBUG(){
        double xmax=0;
        for (CellOL cell : this) {
            if(cell.Xpt()>xmax & cell.type==CellOL.CellType.tumor) {
                xmax = cell.Xpt();
            }
            if(cell.Xpt()>xmax & cell.type==CellOL.CellType.wgd) {
                xmax = cell.Xpt();
                int rvar=0;

            }
        }
        return(xmax);
    }

    public void shiftFront(double x,double y, double xfrac){
        for (CellOL cell : this) {
            if(cell.Xpt()>x*xfrac){
                //CellOL child = this.NewAgentPTSafe(cell.Xpt()-x*(1-xfrac),cell.Ypt());
                CellOL child = this.NewAgentPT(cell.Xpt()-x*(1-xfrac),cell.Ypt());
                child.type=cell.type;
                if(cell.type==CellOL.CellType.vessel) cell.Dispose();
            }else{
            //if(cell.Xpt()<=x*xfrac){
                if(cell.Xpt()>x*(1-xfrac)){
                    cell.MoveSafePT(cell.Xpt()-x*(1-xfrac),cell.Ypt());
                } else{
                    cell.Dispose();
                }
            }


        }

    }

}

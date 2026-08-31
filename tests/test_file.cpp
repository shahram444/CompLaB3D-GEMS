#include <cstdio>
#include <cmath>
#include "complab3d_symbolic.hh"
using namespace complab_sym;
int fails = 0;
void ck(const char *w, double g, double t, double tol=1e-9){
    bool ok = std::fabs(g-t) <= tol*(1.0+std::fabs(t)); if(!ok) ++fails;
    std::printf("  %-44s got %-13.7g want %-13.7g %s\n", w,g,t, ok?"ok":"** FAIL **"); }

int main(){
    Program P; std::string err;
    if(!load(P,"../examples/16_learned_rate_laws/fitting/ecoli.sym",&err)){ std::printf("LOAD FAILED: %s\n",err.c_str()); return 1; }
    std::printf("%s", describe(P,"../examples/16_learned_rate_laws/fitting/ecoli.sym").c_str());

    std::printf("\n--- structure\n");
    std::printf("  %-44s %d\n","variables",(int)P.vars.size());
    std::printf("  %-44s %d\n","rates",(int)P.rates.size());
    std::printf("  %-44s %s\n","rate order", (P.rates[0].name+", "+P.rates[1].name+", "+P.rates[2].name).c_str());

    std::printf("\n--- evaluation, inside the range\n");
    std::vector<double> v; v.push_back(5.0); v.push_back(10.0); v.push_back(1.0);
    std::vector<double> out; long cl=0;
    evaluate(P,v,out,true,&cl);
    double growth = 0.0412*5.0*10.0/(0.083+10.0);
    ck("growth",  out[0], growth);
    ck("acetate = 0.31*growth (chained)", out[1], 0.31*growth);
    ck("glucose uptake (chained)", out[2], -growth/0.0412);
    ck("clamped count", (double)cl, 0);

    std::printf("\n--- clamping outside the fitted box\n");
    std::vector<double> hi; hi.push_back(50.0); hi.push_back(10.0); hi.push_back(1.0);
    cl=0; evaluate(P,hi,out,true,&cl);
    double atTen = 0.0412*10.0*10.0/(0.083+10.0);
    ck("glucose 50 clamps to 10", out[0], atTen);
    ck("clamp counted", (double)cl, 1);
    cl=0; evaluate(P,hi,out,false,&cl);
    ck("clamp off -> extrapolates", out[0], 0.0412*50.0*10.0/(0.083+10.0));
    ck("nothing counted", (double)cl, 0);

    std::printf("\n--- Ecoli has no range line, so it is never clamped\n");
    std::vector<double> big; big.push_back(5.0); big.push_back(10.0); big.push_back(1e6);
    cl=0; evaluate(P,big,out,true,&cl);
    ck("no clamp for an unranged variable",(double)cl,0);

    std::printf("\n--- bad files are refused with a useful message\n");
    const char *bad[][2] = {
      {"b1.sym","vars a\nrate x = y*2\n"},
      {"b2.sym","rate x = 1\n"},
      {"b3.sym","vars a\nrate x = a*\n"},
      {"b4.sym","vars a\nrate x = a\nrate x = 2*a\n"},
      {"b5.sym","vars a\nrange b 1 2\n"},
      {"b6.sym","vars a\nrange a 5 1\n"},
      {"b7.sym","vars a\nwibble 3\n"},
      {"b8.sym","vars a\nrate x  a*2\n"},
      /* forward reference: acetate uses growth before growth exists */
      {"b9.sym","vars a\nrate acetate = 0.3*growth\nrate growth = a\n"} };
    for(int i=0;i<9;++i){
        std::FILE*f=std::fopen(bad[i][0],"w"); std::fputs(bad[i][1],f); std::fclose(f);
        Program Q; std::string e;
        bool okLoad = load(Q,bad[i][0],&e);
        std::printf("  %-10s %s\n", bad[i][0], okLoad? "** FAIL: accepted **" : ("refused: "+e).c_str());
        if(okLoad) ++fails;
    }
    std::printf("\n%s\n", fails? "SOME CHECKS FAILED":"all file checks passed");
    return fails;
}

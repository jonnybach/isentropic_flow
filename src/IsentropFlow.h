#ifndef ISENTROPFLOW_H_
#define ISENTROPFLOW_H_

static const double PI = 3.141592653589;
static const double Rgas = 286.9; //air gas constant J/kg-K

typedef enum {
	kFlowSubsonic = 0,
	kFlowSupersonic = 1
} FlowRegime;

double VATfM(const double gamma, const double Mach);
double VELfTTM(const double gamma, const double Tt, const double Mach);
double QoverPT(const double gamma, const double Mach);
double PMturn_NUfM(const double gamma, const double Mach);
double PMturn_MfNU(const double gamma, const double Nu);
double MfNu(const double gamma, const double Nu);
double MfNu(const double gamma, const double Nu);
double Mu2(const double gamma, const double Mach);
double P2P1_PMexpansion(const double gamma, const double Mach, const double delta_ang_deg);
double TRfM(const double gamma, const double Mach);
double MfTR(const double gamma, const double Tr);
double PRfM(const double gamma, const double Mach);
double MfPR(const double gamma, const double Pr);
double DRfM(const double gamma, const double Mach);
double MfDR(const double gamma, const double Dr);
double FPfM(const double gamma, const double Mach);
double MfFP(const double gamma, const double wtap, const FlowRegime flowSpeed);
double ARfM(const double gamma, const double Mach);
double MfAR(const double gamma, const double Aratio, const FlowRegime flowSpeed);

double FSAOAS(const double gamma, const double Aratio);
double FOPTA(const double gamma, const double NPR);
double FOWT(const double gamma, const double NPR);
double Cv(const double Cs, const double gamma, const double AeOAs, const double NPR);
double CvfCs(const double Cs, const double Cd, const double gamma
        , const double ARatio, const double NPR);
double CsfCv(const double CV, const double Cd, const double gamma
	, const double ARatio, const double NPR);

double CD_ASME(const double ReD);

#endif /*ISENTROPFLOW_H_*/

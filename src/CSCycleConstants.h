/*
 * CSCycleConstants.h
 *
 *  Created on: Sep 19, 2018
 *      Author: dhf
 *      Constants relating to carbonate-silicate cycle implementation
 *      (Haqq-Misra et al 2016)
 */

#ifndef CSCYCLECONSTANTS_H_
#define CSCYCLECONSTANTS_H_


const double CO2Earth = 3.3e-4; // Earth's current CO2 partial pressure
const double kactive = 0.09;
const double krun = 0.045;
const double outgassingRateEarth = 7.0e-8/year; // outgassing rate in bars per second

const double aLand = 0.2;        // surface albedo of land
const double aIceVisible = 0.8;  // surface albedo of water ice in optical
const double aIceIR = 0.5;       // surface albedo of water ice in IR
const double aCO2Ice = 0.35;     // surface albedo of CO2 ice
const double fCloud = 0.07;      // Fraction of H2O cloud cover

// Fraction of spectrum in visible for F,G,K,M stars

const double fStarFVisible = 0.67;
const double gStarFVisible = 0.52;
const double kStarFVisible = 0.32;
const double mStarFVisible = 0.1;

const double TsatSolid = 216.56;

// IR cooling coefficients

const double irco[] =
    {
	    9.12805643734088612007e+00,
	    4.58408794776474781685e+00,
	    -8.47261075511868995136e+01,
	    +4.35517381110739065786e-01,
	    -2.86355036266497364750e+01,
	    +2.96626642450453971378e+02,
	    -6.01082900358798077889e-02,
	    -2.60414691486032312540e+00,
	    +5.69812976578495309354e+01,
	    -4.62596100050751886101e+02,
	    +2.18159373001554097310e-03,
	    +1.61456772400849241089e-01,
	    +3.75623788186383888998e+00,
	    -3.53347289235212116409e+01,
	    +2.75011005363597746509e+02
    };

const vector<double> ircoeff(irco, irco + sizeof(irco) / sizeof(irco[0]));

// Ocean albedo developed from Fresnel tables
// albco[i] = ocean reflectance at i degrees

const double albco[] =
	{
		0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
		0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021, 0.021,
		0.021, 0.021, 0.021, 0.021, 0.021, 0.022, 0.022, 0.022, 0.022, 0.022,
		0.022, 0.023, 0.023, 0.023, 0.023, 0.024, 0.024, 0.024, 0.024, 0.025,
		0.025, 0.025, 0.026, 0.027, 0.028, 0.029, 0.030, 0.031, 0.032, 0.034,
		0.035, 0.036, 0.038, 0.040, 0.042, 0.044, 0.047, 0.050, 0.054, 0.058,
		0.062, 0.066, 0.070, 0.075, 0.082, 0.096, 0.104, 0.114, 0.114, 0.124,
		0.136, 0.148, 0.162, 0.178, 0.196, 0.215, 0.238, 0.260, 0.288, 0.314,
		0.350, 0.386, 0.428, 0.476, 0.529, 0.586, 0.650, 0.720, 0.806, 0.896
	};


const vector<double> oceanAlbedo(albco, albco + sizeof(albco) / sizeof(albco[0]));

// Albedo Coefficients for M Stars

const double mStarHot[] =
    {
	    -7.15889280368259939280e-01,
	    -1.98136800730935841441e-01,
	    +9.84775449827866244945e-01,
	    -1.56391572586811709866e-02,
	    -6.27978666433394616675e-01,
	    +2.72826541972741240527e-02,
	    -1.83968255888854792524e+00,
	    -1.76445542493349726010e-02,
	    +4.88122504773761800578e+00,
	    +8.16466305082373011714e+00,
	    -9.42628034611611292926e-02,
	    -4.09511921617251815064e+01,
	    -4.28369698215551460933e-04,
	    +2.55820387203858357061e-01,
	    +4.98181255672654828004e+01,
	    +1.74804687259385131692e-02,
	    -1.78959455813329754159e-01,
	    +9.20933759574703303397e-03,
	    +4.33889186559404560484e-01,
	    +4.94168389369200689032e+00,
	    +4.28308713193277179609e-01,
	    -2.57914484202699831883e+01,
	    -9.90612736197062256072e-03,
	    -1.12652187217021260146e+00,
	    +3.37025608178292372941e+01,
	    +4.95822115912239240743e+01,
	    +5.86747385748280669837e-01,
	    -3.75850672216186808328e+02,
	    +8.26245532595887566263e-02,
	    -2.52798943533546927043e+00,
	    +9.48518044641436858910e+02,
	    +7.58391275375584609543e-04,
	    -1.94943546896162750448e-01,
	    +2.68713810507689165874e+00,
	    -7.96510574573562053047e+02,
    };

const double mStarCold[] =
    {
	    -5.32746030999423081376e-01,
	    -7.19022284759715302194e-01,
	    +5.03304374139034726987e-01,
	    +1.40117784903560482074e-03,
	    +5.20308913208579193466e-01,
	    +1.84649470604131248075e-02,
	    -4.30909214152957076305e-01,
	    -2.52209479270658273875e-02,
	    +2.04657319578802932014e+00,
	    -2.24771771283140209263e+00,
	    +1.73986782860104903514e-02,
	    +9.97383536921341118386e+00,
	    -4.44911599245020991829e-05,
	    -2.97880877255313897267e-02,
	    -1.26216890838289153010e+01,
	    +2.00320513295134793735e-03,
	    +5.33365308997984510775e-02,
	    +6.02794793065432234908e-03,
	    -1.03172773926871935712e-01,
	    -7.20688958123098100117e+00,
	    +1.84702956139857557560e-02,
	    +3.26515510288737047517e+01,
	    -1.44154978290867378382e-02,
	    -1.45491367231501134150e-01,
	    -3.66323399007161611962e+01,
	    -2.14507917132875682853e+01,
	    +3.47529790408816263714e-01,
	    +1.50422575721380695768e+02,
	    +1.24670173663235066275e-02,
	    -1.52798822618099761073e+00,
	    -3.51237344120818193005e+02,
	    -1.18980466901116377701e-03,
	    -3.52156741076393511869e-02,
	    +1.66607899218190680379e+00,
	    +2.73570694536131782115e+02
    };

const vector<double> mStarHotAlbedo(mStarHot,mStarHot + sizeof(mStarHot) / sizeof(mStarHot[0]));

const vector<double> mStarColdAlbedo(mStarCold, mStarCold + sizeof(mStarCold) / sizeof(mStarCold[0]));

// Albedo coefficients for K Stars

const double kStarHot[] =
    {
	    -5.90615904550030323961e-01,
	    -2.37208462719700086119e-01,
	    +1.13931290171631860453e+00,
	    -3.08377226984211140481e-02,
	    -1.28102871876537105500e+00,
	    +6.19626526274854941279e-02,
	    -1.69774543632547247896e+00,
	    -2.04235309152326696691e-02,
	    +4.57773963051313170780e+00,
	    +4.12579981898720138389e+00,
	    -4.71535899181005566105e-02,
	    -2.11943079255920068249e+01,
	    +1.04512904730470691968e-03,
	    +1.63558516256956165691e-01,
	    +2.58180027522388897410e+01,
	    +6.73502603995534410153e-02,
	    -2.86882220125632747543e-01,
	    +1.86555392184974862257e-02,
	    +6.39804770045043325055e-01,
	    +6.68592076399450507829e-01,
	    +3.40006379075857678718e-01,
	    -4.47130851866186240784e+00,
	    -1.37881811958935548978e-02,
	    -9.37032157741380156146e-01,
	    +7.18430533034758767030e+00,
	    +3.55789690699431560006e+01,
	    +6.05324898846670844677e-01,
	    -2.69163119251510408958e+02,
	    +7.16411109217734848320e-02,
	    -2.70086156378175923365e+00,
	    +6.77928855364022979302e+02,
	    +1.75159795783044118338e-03,
	    -1.57430537384521729294e-01,
	    +3.03311085863799512197e+00,
	    -5.67991643797108054059e+02
    };

const double kStarCold[] =
    {
	    -4.64200702979755630562e-01,
	    -7.08632250700805155219e-01,
	    +2.94280811500949868176e-01,
	    -7.03623530512023729472e-03,
	    +8.32291767736463627969e-01,
	    +4.17366696114838467424e-02,
	    -3.27381031519621223946e-01,
	    -2.01262118534941206183e-02,
	    +1.79927783353684689338e+00,
	    -1.52517039855256109071e+00,
	    +1.37559348494792656192e-02,
	    +6.83260807978973971899e+00,
	    -8.41179377240318719661e-06,
	    -1.49071820222227719582e-02,
	    -9.12614113110200442236e+00,
	    +6.00961538795274658603e-03,
	    +5.37129271957198853316e-02,
	    +1.53045807573081431990e-02,
	    -7.29624602174949454803e-02,
	    -5.33166179727468403371e+00,
	    -1.88071318117267599623e-02,
	    +2.40599196073251704320e+01,
	    -1.45876963778264045341e-02,
	    -6.74760640962892566108e-02,
	    -2.68205138267828857579e+01,
	    -1.50036243805226146009e+01,
	    +3.43803701611855749842e-01,
	    +1.05678253976713648399e+02,
	    +8.82128516804767635884e-03,
	    -1.51838731559315354147e+00,
	    -2.47778019102939651930e+02,
	    -1.97405245617707729823e-04,
	    -1.65222744277499604404e-02,
	    +1.69889556964086874125e+00,
	    +1.93919792517779910668e+02
    };

const vector<double> kStarHotAlbedo(kStarHot, kStarHot + sizeof(kStarHot) / sizeof(kStarHot[0]));
const vector<double> kStarColdAlbedo(kStarCold, kStarCold + sizeof(kStarCold) / sizeof(kStarCold[0]));

const double gStarHot[] =
    {
	    -4.41391619954555503025e-01,
	    -2.60017516002879089942e-01,
	    +1.08110772295329837789e+00,
	    -3.93863285843020910493e-02,
	    -1.46383456258096611435e+00,
	    +9.91383778608142668398e-02,
	    -1.45914724229303338632e+00,
	    -2.72769392852398387395e-02,
	    +3.99933641081463919775e+00,
	    +1.07231336256525633388e+00,
	    -1.04302520934751417891e-02,
	    -6.10296439299006454604e+00,
	    +2.69255203910960137434e-03,
	    +9.50143253373007257157e-02,
	    +7.37864215757422226005e+00,
	    +1.28580729156335171748e-01,
	    -3.07800300913486257759e-01,
	    +2.27715594632176554502e-02,
	    +6.11699085276039222769e-01,
	    -2.33213409642421742873e+00,
	    +2.56011431303802661219e-01,
	    +1.05912148222549546972e+01,
	    -1.85772688884413561539e-02,
	    -7.55796861024326749323e-01,
	    -1.16485004141808623501e+01,
	    +2.74062491988752192640e+01,
	    +5.46044240911252587445e-01,
	    -2.05761674358916081928e+02,
	    +5.57943359123403426203e-02,
	    -2.49880329758542751861e+00,
	    +5.14448995054491206247e+02,
	    +2.43702089287719950508e-03,
	    -1.09384840764980617589e-01,
	    +2.92643187434628071486e+00,
	    -4.27802454850920923946e+02
    };

const double gStarCold[] =

    {
	    -3.64301272050786051349e-01,
	    -6.66571453035937344644e-01,
	    +1.38761634791769922215e-01,
	    -1.40826323888164368220e-02,
	    +9.41440608298288128530e-01,
	    +7.10961643487220129600e-02,
	    -2.19180456421237290776e-01,
	    -1.82873271476295846949e-02,
	    +1.48505536251773073708e+00,
	    -9.01309617860975631487e-01,
	    +1.92113767482554841093e-02,
	    +4.11334031794617160926e+00,
	    +6.80906172782627400891e-04,
	    -1.66632232847024261413e-02,
	    -6.01321219414692986760e+00,
	    +5.20833333338503734478e-02,
	    +1.09511892935421337181e-01,
	    +1.86369741605604787027e-02,
	    -2.54092206932019781807e-01,
	    -4.00290429315177131997e+00,
	    -4.60694421170402754195e-02,
	    +1.79103047870275950970e+01,
	    -1.59834667195196747369e-02,
	    -1.29954198131196525801e-02,
	    -1.97041106668471570629e+01,
	    -9.28987827590191805882e+00,
	    +2.33079221557892068972e-01,
	    +6.58750181054108310263e+01,
	    +7.46763857253681870296e-03,
	    -1.00561681124449076030e+00,
	    -1.55355955538023465579e+02,
	    +7.11268878229609079374e-04,
	    -3.36136500021004319336e-03,
	    +1.13977221457453326003e+00,
	    +1.22439629486842392225e+02,
    };

const vector<double> gStarHotAlbedo(gStarHot, gStarHot + sizeof(gStarHot) / sizeof(gStarHot[0]));
const vector<double> gStarColdAlbedo(gStarCold, gStarCold + sizeof(gStarCold) / sizeof(gStarCold[0]));

const double fStarHot[] =
    {
	    -2.94907643999552049330e-01,
	    -5.86780253976008414618e-01,
	    +6.49675557670263970067e-02,
	    -1.79182590979124686803e-02,
	    +8.94400119847050145694e-01,
	    +1.02431720374415072272e-01,
	    -1.95479829654055359267e-01,
	    -2.12788556357093004701e-02,
	    +1.30481589291653921059e+00,
	    -3.67440225715104984427e-01,
	    +1.79968532967263519784e-02,
	    +1.75307643417233238736e+00,
	    +1.70822171015460829747e-03,
	    -1.20533628790382989389e-03,
	    -3.23742655000626289308e+00,
	    +8.94764957274146316424e-02,
	    +1.06659840879757411569e-01,
	    +1.97277310982091197422e-02,
	    -2.90643338479389279350e-01,
	    -3.02285149426790678007e+00,
	    -6.27300427194927429086e-02,
	    +1.34578175237932615715e+01,
	    -1.78176771900690443518e-02,
	    +1.66678696068688798892e-02,
	    -1.46647101613171422230e+01,
	    -4.01170875081073852186e+00,
	    +1.27577468774495134118e-01,
	    +2.89809987786521858766e+01,
	    +5.53139586753823280646e-03,
	    -5.22907047273701075518e-01,
	    -6.94592466763898386262e+01,
	    +1.50903897147278446764e-03,
	    +9.57048015146020475408e-03,
	    +6.14315422188238868806e-01,
	    +5.58749555497795213910e+01,
    };

const double fStarCold[] =
    {
    -3.31402546093767180757e-01,
    -2.47134807241273235512e-01,
    +8.81397154521980197295e-01,
    -4.16225627219204072360e-02,
    -1.23908968458793888878e+00,
    +1.30930881442010454974e-01,
    -1.20262852602715497596e+00,
    -3.59480490797699711591e-02,
    +3.32522395760313171920e+00,
    -4.71338194513966635135e-01,
    +1.09104030049774064492e-02,
    +1.69353435116822814699e+00,
    +4.39360682804632116799e-03,
    +5.93274962992209262125e-02,
    -2.26669461112944992109e+00,
    +1.68619791658343592955e-01,
    -2.99057302813440284428e-01,
    +2.12529784002843634416e-02,
    +5.38685716278799464618e-01,
    -4.01162274610886449011e+00,
    +1.98990869102443984628e-01,
    +1.91535828079575587424e+01,
    -2.20977011730780145693e-02,
    -6.28983639418757678463e-01,
    -2.25212242581750210491e+01,
    +1.80081287011609418869e+01,
    +4.64223667003067608228e-01,
    -1.34086408527571904870e+02,
    +4.14114309723633972071e-02,
    -2.16873758498834590114e+00,
    +3.32495385798052097925e+02,
    +2.75346362673409110478e-03,
    -6.88053255006287611817e-02,
    +2.62954623889728411612e+00,
    -2.73995632394573760848e+02
    };

const vector<double> fStarHotAlbedo(fStarHot, fStarHot+sizeof(fStarHot)/sizeof(fStarHot[0]));
const vector<double> fStarColdAlbedo(fStarCold, fStarCold+sizeof(fStarCold)/sizeof(fStarCold[0]));


#endif /* CONSTANTS_H_ */
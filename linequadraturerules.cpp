#include "linequadraturerules.h"

lineQuadratureRules::lineQuadratureRules()
{

}

Eigen::MatrixX2d lineQuadratureRules::weightsandAbscissa(int n) // for integrating over [-1, 1]
{
    Eigen::MatrixX2d abscissaAndWeights;
    if(n <= 0)
    {
        n = 1;
    }
    if(n == 1)
    {
        abscissaAndWeights = Eigen::MatrixX2d(1,2);

        abscissaAndWeights.col(0) <<
        0.0;

        abscissaAndWeights.col(1) <<
        2.0;
    }
    if(n == 2)
    {
        abscissaAndWeights = Eigen::MatrixX2d(2,2);

        abscissaAndWeights.col(0) <<
        -0.577350269189625764509148780502,
        0.577350269189625764509148780502;

        abscissaAndWeights.col(1) <<
        1.0000000000000000,
        1.0000000000000000;
    }
    if(n == 3)
    {
        abscissaAndWeights = Eigen::MatrixX2d(3,2);

        abscissaAndWeights.col(0) <<
        -0.774596669241483377035853079956,
        0.000000000000000000000000000000,
        0.774596669241483377035853079956;

        abscissaAndWeights.col(1) <<
         0.5555555555555556,
         0.8888888888888889,
         0.5555555555555556;
    }
    if(n == 4)
    {
        abscissaAndWeights = Eigen::MatrixX2d(4,2);

        abscissaAndWeights.col(0) <<
        -0.861136311594052575223946488893,
        -0.339981043584856264802665759103,
        0.339981043584856264802665759103,
        0.861136311594052575223946488893;

        abscissaAndWeights.col(1) <<
        0.347854845137453857373063949222,
        0.652145154862546142626936050778,
        0.652145154862546142626936050778,
        0.347854845137453857373063949222;
    }
    if(n == 5)
    {
        abscissaAndWeights = Eigen::MatrixX2d(5,2);

        abscissaAndWeights.col(0) <<
        -0.906179845938663992797626878299,
        -0.538469310105683091036314420700,
         0.000000000000000000000000000000,
         0.538469310105683091036314420700,
         0.906179845938663992797626878299;

        abscissaAndWeights.col(1) <<
        0.236926885056189087514264040720,
        0.478628670499366468041291514836,
        0.568888888888888888888888888889,
        0.478628670499366468041291514836,
        0.236926885056189087514264040720;
    }
    if(n > 5)
    {
        Eigen::VectorXd xAbszissa(16);
        xAbszissa << 0.0483076656877383162348126,0.1444719615827964934851864,0.2392873622521370745446032,0.3318686022821276497799168,0.4213512761306353453641194,0.5068999089322293900237475,0.5877157572407623290407455,0.6630442669302152009751152,0.7321821187402896803874267,0.7944837959679424069630973,0.8493676137325699701336930,0.8963211557660521239653072,0.9349060759377396891709191,0.9647622555875064307738119,0.9856115115452683354001750,0.9972638618494815635449811;
        Eigen::VectorXd xWeights(16);
        xWeights << 0.0965400885147278005667648,0.0956387200792748594190820,0.0938443990808045656391802,0.0911738786957638847128686,0.0876520930044038111427715,0.0833119242269467552221991,0.0781938957870703064717409,0.0723457941088485062253994,0.0658222227763618468376501,0.0586840934785355471452836,0.0509980592623761761961632,0.0428358980222266806568786,0.0342738629130214331026877,0.0253920653092620594557526,0.0162743947309056706051706,0.0070186100094700966004071;

        abscissaAndWeights = Eigen::MatrixX2d(32,2);
        abscissaAndWeights.col(0).head(16) = -xAbszissa;
        abscissaAndWeights.col(0).tail(16) = xAbszissa;
        abscissaAndWeights.col(1).head(16) = xWeights;
        abscissaAndWeights.col(1).tail(16) = xWeights;
    }
    return abscissaAndWeights;
}
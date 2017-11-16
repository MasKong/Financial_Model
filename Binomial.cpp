#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//# Author: Mingkang LI
//# Binomial tree (Cox-Rox-Rubenstein) for European/American Option Valuation
//
//# Inspired on:
//# Binomial Tree for America and European options by Mehdi Bounouar
//# Binomial Tree Option Valuation Cox, Ross, Rubinstein method by "www.quantandfinancial.com"

class Binomial{
public:
    double S,X,T,sigma,r;           //S is current stock price, X is option exercise price
    unsigned long n;                          //T is option life expiration time in year, sigma is volatility
    char trade_type;                  //r is risk free interest rate, n is step in binomial tree


protected:
    double p,u,d,delta_t;
    vector<double> stock_value, option_value;
    vector<double> forward();
//    virtual vector<double> backward();
//    vector<double> backward();
//    vector<double> backward_put();
//    double calculate();
    double max_E(double v1, double v2), max_A(double v1, double v2, double v3);

};

class European : Binomial{
public:
//    European();

    European(double SS, double XX,double TT,
             double sigmas,double rr, unsigned long nn, char trade);
//    double calculate();


private:
    vector<double> backward();
    vector<double> backward_put();
};

class American : Binomial{
public:

    American(double SS, double XX,double TT,
             double sigmas,double rr, unsigned long nn, char trade);
    double calculate();


private:
    vector<double> backward();
    vector<double> backward_put();
};

European::European(double SS, double XX,double TT,
                   double sigmas,double rr, unsigned long nn, char trade)
//:S{S},X{X},T{T},sigma{sigma}, r{r}, n{n}
{
//    p = exp()
    S = SS;
    X = XX;
    T = TT;
    sigma = sigmas;
    r = rr;
    n = nn;
    trade_type = trade;
    delta_t = T/n;
    u = exp(sigma * sqrt(delta_t));
    d = 1/u;
    p = (exp(r*delta_t)-d)/(u-d);

}

American::American(double SS, double XX,double TT,
                   double sigmas,double rr, unsigned long nn, char trade)
//:S{S},X{X},T{T},sigma{sigma}, r{r}, n{n}
{
//    p = exp()
    S = SS;
    X = XX;
    T = TT;
    sigma = sigmas;
    r = rr;
    n = nn;
    trade_type = trade;
    delta_t = T/n;
    u = exp(sigma * sqrt(delta_t));
    d = 1/u;
    p = (exp(r*delta_t)-d)/(u-d);

}

double Binomial::max_E(double v1, double v2) {
    return (v1>v2)?v1:v2;

}

double Binomial::max_A(double v1, double v2, double v3) {
    double max_2 = (v1>v2)?v1:v2;
    return (v3>max_2)?v3:max_2;

}


vector<double> Binomial::forward(){        //compute the stock value for each step
    unsigned long total = (n+1)*(n+2)/2;      //total number of stock value and option value
    vector<double> stock_value(total);
    stock_value[0] = S;
    int step = 1;
    long value_num = (step+1)*(step+2)/2;
    for (int i=1; i<total; ++i){
        if (i+1 > (step+1)*(step+2)/2){
            step += 1;
            value_num = (step+1)*(step+2)/2;
        }
        if (i == value_num -1){             //the stock price go down for each step
            stock_value[i] = stock_value[i-step-1] * d;
        }
        else {
            stock_value[i] = stock_value[i - step] * u;
        }
//        cout << stock_value[i] << endl;

    }
    return stock_value;


}


vector<double> European::backward(){
    unsigned long total = (n+1)*(n+2)/2;      //total number of stock value and option value
    vector<double> option_value(total);
    long step = n;

    for (long i=total-1; i>=total-n-1; --i){         //calculate the option value at expiry
        option_value[i] = max_E(stock_value[i]-X,0);


    }
    step -= 1;
    long value_num = (step+0)*(step+1)/2;
    for (long i=total-1-(n+1); i>=0; --i){    //compute backward. There are n+1 elements in n step
        if (i == value_num - 1){
            step -= 1;
            value_num = (step+0)*(step+1)/2;
        }
        option_value[i] = max_E(exp(-r * delta_t) * (p * option_value[i+step+1] + (1-p) * option_value[i+step+2]),0);
//        option_value[i] = exp(-r * delta_t) * (p * option_value[i+step+1] + (1-p) * option_value[i+step+2]);
//        cout << option_value[i] << endl;
    }
    return option_value;
}

vector<double> European::backward_put(){
    unsigned long total = (n+1)*(n+2)/2;      //total number of stock value and option value
    vector<double> option_value(total);
    long step = n;

    for (long i=total-1; i>=total-n-1; --i){         //calculate the option value at expiry
        option_value[i] = max_E(X-stock_value[i],0);


    }
    step -= 1;
    long value_num = (step+0)*(step+1)/2;
    for (long i=total-1-(n+1); i>=0; --i){    //compute backward. There are n+1 elements in n step
        if (i == value_num - 1){
            step -= 1;
            value_num = (step+0)*(step+1)/2;
        }
        option_value[i] = max_E(exp(-r * delta_t) * (p * option_value[i+step+1] + (1-p) * option_value[i+step+2]),0);
//        option_value[i] = exp(-r * delta_t) * (p * option_value[i+step+1] + (1-p) * option_value[i+step+2]);
//        cout << option_value[i] << endl;
    }
    return option_value;
}

//double Binomial::calculate(){
//    stock_value = forward();
//    if (trade_type == 'C')
//        option_value = backward();
//    else
//        option_value = backward_put();
//    return option_value[0];
//}

double American::calculate(){
    stock_value = forward();
    if (trade_type == 'C')
        option_value = backward();
    else
        option_value = backward_put();
    return option_value[0];
}

vector<double> American::backward(){
    unsigned long total = (n+1)*(n+2)/2;      //total number of stock value and option value
    vector<double> option_value(total);
    long step = n;

    for (long i=total-1; i>=total-n-1; --i){         //calculate the option value at expiry
        option_value[i] = max_E(stock_value[i]-X,0);


    }
    step -= 1;
    long value_num = (step+0)*(step+1)/2;
    for (long i=total-1-(n+1); i>=0; --i){    //compute backward. There are n+1 elements in n step
        if (i == value_num - 1){
            step -= 1;
            value_num = (step+0)*(step+1)/2;
        }
        option_value[i] = max_A( exp(-r * delta_t) * (p * option_value[i+step+1] + (1-p) * option_value[i+step+2]),
                                 X-stock_value[i], 0);
//        cout << option_value[i] << endl;
    }
    return option_value;
}

vector<double> American::backward_put(){
    unsigned long total = (n+1)*(n+2)/2;      //total number of stock value and option value
    vector<double> option_value(total);
    long step = n;

    for (long i=total-1; i>=total-n-1; --i){         //calculate the option value at expiry
        option_value[i] = max_E(X-stock_value[i],0);


    }
    step -= 1;
    long value_num = (step+0)*(step+1)/2;
    for (long i=total-1-(n+1); i>=0; --i){    //compute backward. There are n+1 elements in n step
        if (i == value_num - 1){
            step -= 1;
            value_num = (step+0)*(step+1)/2;
        }
        option_value[i] = max_A( exp(-r * delta_t) * (p * option_value[i+step+1] + (1-p) * option_value[i+step+2]),
                                 X-stock_value[i], 0);
//        option_value[i] = exp(-r * delta_t) * (p * option_value[i+step+1] + (1-p) * option_value[i+step+2]);
//        cout << option_value[i] << endl;
    }
    return option_value;
}


int main(){
//    European et(70, 80, 0.25, 0.2, 0.12, 3,'P');
    American at(70, 80, 0.25, 0.2, 0.12, 3,'P');
    double price = at.calculate();
    cout << price << endl;
    return 0;
}




#ifndef WELFORD_H
#define WELFORD_H

class Welford_Algo
{
private:
    double sum_diff_cur;
    double sum_diff_prev=0.0;

    double sum;
    double sum_prev=0.0;
    int num=0;

    void init(double in)
    {
        sum = in;
        ++num;
    }
public:
    Welford_Algo(){}

    void push(double in)
    {
        if(num == 0)
        {
            init(in);
        }
        else
        {
            sum_prev = sum;
            sum += in;
            ++num;

            sum_diff_cur = sum_diff_prev + ( in-sum_prev/(num-1) ) * ( in-sum/num );
            sum_diff_prev = sum_diff_cur;

        }
    }

    double getAverage()
    {
        return sum / num;
    }

    double getVariance()
    {
        return sum_diff_cur/num;
    }
};

#endif // WELFORD_H

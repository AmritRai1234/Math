//! Financial mathematics: interest, NPV, IRR, amortization.

use crate::powers;
use crate::arithmetic::abs;

/// Simple interest: P * r * t.
pub fn simple_interest(principal:f64,rate:f64,time:f64)->f64{principal*rate*time}

/// Compound interest: P * (1+r/n)^(n*t).
pub fn compound_interest(principal:f64,rate:f64,n:f64,time:f64)->f64{
    principal*powers::pow(1.0+rate/n,n*time)
}

/// Continuous compounding: P * e^(r*t).
pub fn continuous_compounding(principal:f64,rate:f64,time:f64)->f64{
    principal*crate::exponential::exp(rate*time)
}

/// Future value of annuity: PMT * ((1+r)^n - 1) / r.
pub fn future_value_annuity(pmt:f64,rate:f64,periods:f64)->f64{
    pmt*(powers::pow(1.0+rate,periods)-1.0)/rate
}

/// Present value of annuity: PMT * (1 - (1+r)^(-n)) / r.
pub fn present_value_annuity(pmt:f64,rate:f64,periods:f64)->f64{
    pmt*(1.0-powers::pow(1.0+rate,-periods))/rate
}

/// Net Present Value: Σ CF_t / (1+r)^t.
pub fn npv(rate:f64,cashflows:&[f64])->f64{
    cashflows.iter().enumerate().map(|(t,&cf)|cf/powers::pow(1.0+rate,t as f64)).sum()
}

/// Internal Rate of Return (Newton's method).
pub fn irr(cashflows:&[f64])->f64{
    let mut r=0.1;
    for _ in 0..1000{
        let pv:f64=cashflows.iter().enumerate().map(|(t,&cf)|cf/powers::pow(1.0+r,t as f64)).sum();
        let dpv:f64=cashflows.iter().enumerate().map(|(t,&cf)|-(t as f64)*cf/powers::pow(1.0+r,t as f64+1.0)).sum();
        if abs(dpv)<1e-30{break;}
        let new_r=r-pv/dpv;
        if abs(new_r-r)<1e-12{return new_r;}
        r=new_r;
    }
    r
}

/// Monthly mortgage payment.
pub fn mortgage_payment(principal:f64,annual_rate:f64,years:u32)->f64{
    let r=annual_rate/12.0;let n=(years*12)as f64;
    principal*r*powers::pow(1.0+r,n)/(powers::pow(1.0+r,n)-1.0)
}

/// Amortization schedule: returns Vec of (payment, principal, interest, balance).
pub fn amortization_schedule(principal:f64,annual_rate:f64,years:u32)->Vec<(f64,f64,f64,f64)>{
    let pmt=mortgage_payment(principal,annual_rate,years);
    let r=annual_rate/12.0;let n=years*12;
    let mut balance=principal;let mut schedule=Vec::with_capacity(n as usize);
    for _ in 0..n{let interest=balance*r;let princ=pmt-interest;balance-=princ;
        schedule.push((pmt,princ,interest,balance.max(0.0)));}
    schedule
}

/// Rule of 72: approximate years to double at given rate.
pub fn rule_of_72(rate_percent:f64)->f64{72.0/rate_percent}

/// CAGR: Compound Annual Growth Rate.
pub fn cagr(begin_value:f64,end_value:f64,years:f64)->f64{
    powers::pow(end_value/begin_value,1.0/years)-1.0
}

/// Payback period (simple): years until cumulative cashflows >= investment.
pub fn payback_period(investment:f64,annual_cashflow:f64)->f64{investment/annual_cashflow}

/// Discounted payback period.
pub fn discounted_payback(investment:f64,cashflows:&[f64],rate:f64)->f64{
    let mut cum=0.0;
    for (t,&cf) in cashflows.iter().enumerate(){
        cum+=cf/powers::pow(1.0+rate,t as f64+1.0);
        if cum>=investment{return(t+1)as f64;}
    }
    f64::INFINITY
}

/// Present value.
pub fn present_value(future_value:f64,rate:f64,periods:f64)->f64{
    future_value/powers::pow(1.0+rate,periods)
}

/// Future value.
pub fn future_value(present_val:f64,rate:f64,periods:f64)->f64{
    present_val*powers::pow(1.0+rate,periods)
}

#[cfg(test)]
mod tests{use super::*;use crate::arithmetic::approx_eq;
    #[test]fn test_compound(){let fv=compound_interest(1000.0,0.05,12.0,10.0);assert!(fv>1647.0&&fv<1649.0);}
    #[test]fn test_npv(){let cf=[-1000.0,300.0,400.0,500.0,200.0];assert!(npv(0.1,&cf)>100.0);}
    #[test]fn test_irr(){let cf=[-1000.0,300.0,400.0,500.0,200.0];let r=irr(&cf);assert!(r>0.1&&r<0.2);}
    #[test]fn test_mortgage(){let pmt=mortgage_payment(200000.0,0.05,30);assert!(pmt>1070.0&&pmt<1075.0);}
    #[test]fn test_rule72(){assert!(approx_eq(rule_of_72(6.0),12.0,1e-10));}
}

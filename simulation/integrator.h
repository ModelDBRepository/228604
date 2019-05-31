#ifndef __INTEGRATOR_H
#define __INTEGRATOR_H

#include "autoparams.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <stack>
namespace integrator 
{

double infty = 1000;

class WorkspaceHandle;

class Workspaces
{
public:
    static const size_t WORKSPACE_SIZE = 1000000;
    static const size_t NUM_INITIAL_WORKSPACES = 4;

    static WorkspaceHandle get_workspace_handle();
    
    virtual ~Workspaces()
    {
        while (!workspace_stack_.empty()) 
        {
            gsl_integration_workspace_free(workspace_stack_.top());
            workspace_stack_.pop();
        }
    }
private:
    Workspaces(): // private constructor
        workspace_stack_()
    {
        for (int i = 0; i < NUM_INITIAL_WORKSPACES; i++)
            workspace_stack_.push(
                    gsl_integration_workspace_alloc(WORKSPACE_SIZE));
    }
    Workspaces(Workspaces const& copy);            // not impl
    Workspaces& operator=(Workspaces const& copy); // not impl
    std::stack<gsl_integration_workspace*> workspace_stack_;
};

class WorkspaceHandle
{
public:
    WorkspaceHandle(std::stack<gsl_integration_workspace*>* workspaces):
        workspaces_(workspaces), ws_(0)
    {
        if (!workspaces->empty())
        {
            ws_ = workspaces_->top();
            workspaces_->pop();
        }
        else ws_ = gsl_integration_workspace_alloc(Workspaces::WORKSPACE_SIZE);
        
    }

    virtual ~WorkspaceHandle()
    {
        workspaces_->push(ws_);
    }

    gsl_integration_workspace* get()
    {
        return ws_;
    }
private:
    std::stack<gsl_integration_workspace*>* workspaces_;
    gsl_integration_workspace* ws_;
};

WorkspaceHandle Workspaces::get_workspace_handle()
{
        static Workspaces workspaces;
        return WorkspaceHandle(&workspaces.workspace_stack_);
}    


struct WrapperInfo
{
    double l;
    double r;
    double ts;
    gsl_function fu;
};

double wrapper_func(const double x, const void* p)
{
    const WrapperInfo* wi = static_cast<const WrapperInfo*>(p);
    const double l = wi->l;
    const double r = wi->r;
    const double ts = wi->ts;
    return (r-l)/2 * 1./atan(ts) * 1./(1.+x*x) * 
        wi->fu.function((r-l)/2*atan(x)/atan(ts)+(l+r)/2, wi->fu.params);
}


typedef double (*func_with_nonconst_args)(double, void*);

double integrate(double (*func)(const double, const void*), const double a, const double b,
        const double eps, const void* const p=0) 
{
    WorkspaceHandle w = Workspaces::get_workspace_handle();
    double result = 0;
    double error = 0;
    gsl_function fu = {(func_with_nonconst_args)func, const_cast<void*>(p)};
    gsl_integration_qags(&fu, a, b, params::integ_epsabs, eps, Workspaces::WORKSPACE_SIZE, 
        w.get(), &result, &error);
    return result;
}

double integrate_inf(double (*func)(double, const void*),
        const double eps, const void* const p=0)
{
    WorkspaceHandle w = Workspaces::get_workspace_handle();
    double result = 0;
    double error = 0;
    gsl_function fu = {(func_with_nonconst_args)func, const_cast<void*>(p)};
    gsl_integration_qagi(&fu, params::integ_epsabs, eps, Workspaces::WORKSPACE_SIZE, 
            w.get(), &result, &error);
    return result;
}

double integrate_infu(double (*func)(double, const void*), const double a, 
        const double eps, const void* const p=0)
{
    WorkspaceHandle w = Workspaces::get_workspace_handle();
    double result = 0;
    double error = 0;
    gsl_function fu = {(func_with_nonconst_args)func, const_cast<void*>(p)};
    gsl_integration_qagiu(&fu, a, params::integ_epsabs, eps, Workspaces::WORKSPACE_SIZE, 
            w.get(), &result, &error);
    return result;
}

double integrate_infl(double (*func)(double, const void*), const double b,
        const double eps, const void* const p=0)
{
    WorkspaceHandle w = Workspaces::get_workspace_handle();
    double result = 0;
    double error = 0;
    gsl_function fu = {(func_with_nonconst_args)func, const_cast<void*>(p)};
    gsl_integration_qagil(&fu, b, params::integ_epsabs, eps, Workspaces::WORKSPACE_SIZE, 
            w.get(), &result, &error);
    return result;
}


double edge_emph_integrate(double (*func)(double, const void*), const double l, 
        const double r, const double ts, const double a, const double b, 
        const double eps, const void* const p=0)
{
    WorkspaceHandle w = Workspaces::get_workspace_handle();
    double result = 0;
    double error = 0;
    gsl_function fu = {(func_with_nonconst_args)func, const_cast<void*>(p)};
    WrapperInfo wi = {l, r, ts, fu};
    gsl_function wra_fu = {(func_with_nonconst_args)wrapper_func, &wi};
    gsl_integration_qags(&wra_fu, 
            tan(2./(r-l)*(a-(l+r)/2)*atan(ts)), 
            tan(2./(r-l)*(b-(l+r)/2)*atan(ts)), 
            params::integ_epsabs, eps, Workspaces::WORKSPACE_SIZE,
            //GSL_INTEG_GAUSS21, 
            w.get(), &result, &error);


    return result;
}


};

#endif // __INTEGRATOR_H

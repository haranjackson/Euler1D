double psiSB(double r)
{
    if (r < 0)
        return 0;

    else if (r < 0.5)
        return 2*r;

    else if (r < 1)
        return 1;

    else if (r < 2)
        return r;

    else
        return 2;
}


double psiVL(double r)
{
    return (r<0) ? 0 : 2*r/(1+r);
}


double psiVA(double r)
{
    return (r<0) ? 0 : r*(1+r)/(1+r*r);
}


double psiMB(double r)
{
    if (r < 0)
        return 0;

    else if (r < 1)
        return r;

    else
        return 1;
}

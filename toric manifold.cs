using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Toricmanifold
{

    private FixAngle _alpha;
    private FixAngle _theta;
    private FixAngle _phi;
    private GameObject _target1;
    private GameObject _target2;
    private Camera _main;

    private Vector3 A;
    private Vector3 B;
    private Vector3 AB;
    private float Sx;
    private float Sy;

    //desired screen position 
    private Vector2 pA;
    private Vector2 pB;

    public Toricmanifold(float alpha, float theta, float phi, GameObject target1, GameObject target2)
    {
        _target1 = target1;
        _target2 = target2;
        A = _target1.transform.position;
        B = _target2.transform.position;
        AB = A - B;

        _alpha = new FixAngle(alpha);
        _phi = new FixAngle(phi);
        _theta = new FixAngle(theta);
        //Better way to get the camera?
        _main = Camera.main;
        computeScale();
    }



    public Toricmanifold(float alpha, GameObject target1, GameObject target2) : this(alpha, 0, 0, target1, target2)
    {
    }

    public Vector3 ToWorldPosition()
    {
        Vector3 C = new Vector3(0, 0, 0);
        float last = Mathf.Sin(_alpha.toRad() + _theta.toRad() / 2);


        Vector3 n = -AB;

        //n.projectZ
        Vector2 n2 = new Vector2(n.x, n.y);

        //n2.rotate90();
        float tmp = n2[0];
        n2[0] = -n2[1];
        n2[1] = tmp;

        Vector3 z = new Vector3(n2.x, n2.y, 0);
        Vector3.Normalize(z);
        Vector3 t = Vector3.Cross(z, n);

        Quaternion qT = Quaternion.AngleAxis(_theta.angle(), t);
        Quaternion qP = Quaternion.AngleAxis(_phi.angle(), AB);
        C = A + (qP * qT * AB) * last;
        return C;
    }

    public Quaternion ComputeOrientation(Vector3 campPos, float TiltAngle = 0)
    {
        Quaternion qLook = computeLookAt(campPos);
        Quaternion qTrans = computeLookAtTransition();
        //Tilt works
        Quaternion qPhi = computeTiltAngle(TiltAngle);

        return qPhi * qLook * (Quaternion.Inverse(qTrans));
    }

    private Quaternion computeLookAtTransition()
    {

        Vector2 pO = new Vector2(0, 0);
        Vector2 pM = new Vector2((pA.x+pB.x)/2 , (pA.y + pB.y) / 2);
        Vector3 p3O = new Vector3(0, 0, 1);
        Vector3.Normalize(p3O);
        
        //FIX correct way to calculate Vector3?
        Vector3 p3M = new Vector3(pM.x / Sx, pM.y / Sy, 1);
        Vector3.Normalize(p3M);
        Vector3 m = Vector3.Cross(p3M, p3O);
        float angle = Vector3.Angle(p3M, p3O);

        return Quaternion.AngleAxis(angle, m);
    }

    private Quaternion computeTiltAngle(float tilt)
    {
        FixAngle phi = new FixAngle(tilt);
        return Quaternion.AngleAxis(phi.angle(), _main.transform.forward);
    }


    private Quaternion computeLookAt(Vector3 camPos)
    {
        Vector3 dA = A - camPos;
        Vector3 dB = B - camPos;

        //Enough for ||CA||?
        Vector3 dA2 = dA / dA.magnitude;
        Vector3 dB2 = dB / dB.magnitude;

        //Vector3.Normalize(dA);
        //Vector3.Normalize(dB);

        Vector3 l = 0.5f * (dA2 + dB2);
        return Quaternion.LookRotation(l);
       
    }

    private void computeScale()
    {
        float VerticalfovAngleRAD = _main.fieldOfView * Mathf.Deg2Rad;
        //https://gist.github.com/coastwise/5951291
        float HorizontalfovAngleRAD = 2 *  Mathf.Atan(Mathf.Tan(VerticalfovAngleRAD/2) * _main.aspect);
        Sx = 1 / Mathf.Tan(HorizontalfovAngleRAD / 2);
        Sy = 1 / Mathf.Tan(VerticalfovAngleRAD / 2);
    }

    public void setDesiredPosition(Vector2 screenPos1, Vector2 screenPos2)
    {
        pA = screenPos1;
        pB = screenPos2;
    }
}

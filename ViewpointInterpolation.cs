using System;
using UnityEngine;

public class ViewpointInterpolation
{
    private Toricmanifold _keypoint1;
    private Toricmanifold _keypoint2;
    private float _time0;
    private float _time1;


    public ViewpointInterpolation(Toricmanifold t1, Toricmanifold t2, float time0, float time1)
    {
        _keypoint1 = t1;
        _keypoint2 = t2;
        _time0 = time0;
        _time1 = time1;
    }

    
    public Vector3 finalPosition(float time)
    {
        float gP = findXOverTimePosition(_time0, _time1, time);



        //main paths
        Vector3 firstTrajectory = calculateTrajectory(gP, _keypoint1);
        Vector3 secondTrajectory = calculateTrajectory(gP, _keypoint2);

        //controlling motion along time
        Vector3 position = firstTrajectory * (1 - gP) + secondTrajectory * gP;

        return position;
    }

    public Quaternion finalOrientation(float time, Vector3 position)
    {
        float gF = findXOverTimeOrientation(_time0, _time1, time);

        position = ToricComputing.FromWorldPosition(position,_keypoint1._target1.transform.position, _keypoint1._target2.transform.position);

        Toricmanifold T1 = new Toricmanifold(position.x * Mathf.Rad2Deg, position.y * Mathf.Rad2Deg, position.z * Mathf.Rad2Deg, _keypoint1._target1,  _keypoint1._target2);
        Vector2[] desPosT1= _keypoint1.getDesiredScreenPositions();
        T1.SetDesiredPosition(desPosT1[0],desPosT1[1]);

        position = ToricComputing.FromWorldPosition(position, _keypoint2._target1.transform.position, _keypoint2._target2.transform.position);

        Toricmanifold T2 = new Toricmanifold(position.x * Mathf.Rad2Deg, position.y * Mathf.Rad2Deg, position.z * Mathf.Rad2Deg, _keypoint2._target1, _keypoint2._target2);
        Vector2[] desPosT2 = _keypoint2.getDesiredScreenPositions();
        T2.SetDesiredPosition(desPosT2[0], desPosT2[1]);

        Quaternion qT1 = T1.ComputeOrientation();
        Quaternion qT2 = T2.ComputeOrientation();


        //controlling motion along time
        Quaternion rotation =  Quaternion.Slerp(qT1 ,qT2 , gF);

        return rotation;
    }

    private Vector3 calculateTrajectory(float x, Toricmanifold targetAB)
    {

        Vector3[] targets = { targetAB._target1.transform.position, targetAB._target2.transform.position };
        Vector3 AB = targets[1] - targets[0];

        Vector3 p0 = ToricComputing.FromWorldPosition(_keypoint1.ToWorldPosition(), targets[0], targets[1]);
        Vector3 p1 = ToricComputing.FromWorldPosition(_keypoint2.ToWorldPosition(), targets[0], targets[1]);

     

        float alphaInterpolated = p0.x * (1 - x) + p1.x * x;
        alphaInterpolated *= Mathf.Rad2Deg;

        Vector3 vantageA0 = _keypoint1.ToWorldPosition() - targets[0];
        float distanceA0 = vantageA0.magnitude;
        vantageA0 = vantageA0.normalized;

        Vector3 vantageA1 = _keypoint2.ToWorldPosition() - targets[0];
        float distanceA1 = vantageA1.magnitude;
        vantageA1 = vantageA1.normalized;
        Vector3 interpolatedVa = (1 - x)* vantageA0 + x * vantageA1;
        float interpolatedDistanceA =  (1 - x)* distanceA0 + x * distanceA1;

        Vector3 vantageB0 = _keypoint1.ToWorldPosition() - targets[1];
        float distanceB0 = vantageB0.magnitude;
        vantageB0 = vantageB0.normalized;

        Vector3 vantageB1 = _keypoint2.ToWorldPosition() - targets[1];
        float distanceB1 = vantageB1.magnitude;
        vantageB1 = vantageB1.normalized;
        Vector3 interpolatedVb =  (1 - x)* vantageB0 + x * vantageB1;
        float interpolatedDistanceB =  (1 - x)* distanceB0 + x * distanceB1;

      

        float thetaInterpolatedA = 2 * Vector3.Angle(interpolatedVa, AB);
        float thetaInterpolatedB = 2 * (180 - Vector3.Angle(-AB,interpolatedVb) - alphaInterpolated);
        


        Vector3 upOnPlane = Vector3.ProjectOnPlane(Vector3.up, AB);
        Vector3 VAonPlane = Vector3.ProjectOnPlane(interpolatedVa, AB);
        Vector3 phiZero = Vector3.Cross(upOnPlane, AB);


        float phiInterpolatedA = Vector3.SignedAngle(phiZero, VAonPlane, -AB);


        Toricmanifold intersectionTa = new Toricmanifold( alphaInterpolated , thetaInterpolatedA, phiInterpolatedA, targetAB._target1, targetAB._target2);
       
        float distanceAT = (intersectionTa.ToWorldPosition() - targets[0]).magnitude;


         upOnPlane = Vector3.ProjectOnPlane(Vector3.up, AB);
        Vector3 VBonPlane = Vector3.ProjectOnPlane(interpolatedVa, AB);
         phiZero = Vector3.Cross(upOnPlane, AB);


        float phiInterpolatedB = Vector3.SignedAngle(phiZero, VBonPlane, -AB);

        Toricmanifold intersectionTb = new Toricmanifold(alphaInterpolated, thetaInterpolatedB , phiInterpolatedB, targetAB._target1, targetAB._target2);
       
        float distanceBT = (intersectionTb.ToWorldPosition() - targets[1]).magnitude;


        float lambdaA = Mathf.Sin(Vector3.Angle(AB, interpolatedVa));
        float lambdaB = Mathf.Sin(Vector3.Angle(AB, interpolatedVb));

        Vector3 interpolatedTargets1 = targets[0] + interpolatedVa * 0.5f*(interpolatedDistanceA + distanceAT);
        Vector3 interpolatedeTargets2 = targets[1] + interpolatedVb * 0.5f* (interpolatedDistanceB + distanceBT);
        return 0.5f * ( interpolatedTargets1 + interpolatedeTargets2 );
    }

    private float findXOverTimePosition(float time0, float time1, float time)
    {
        return (time - time0) / (time1-time0);
    }

    private float findXOverTimeOrientation(float time0, float time1, float time)
    {
        return (time - time0) / (time1 - time0);
    }
}

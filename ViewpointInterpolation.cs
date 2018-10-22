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
        keypoint1 = t1;
        keypoint2 = t2;
        _time0 = time0;
        _time1 = time1;
    }


    public Vector3 finalPosition(float time)
    {
        float gP = findXOverTime(_time0, _time1, time);



        //main paths
        Vector3 firstTrajectory = calculateTrajectory(gP, _keypoint1);
        Vector3 secondTrajectory = calculateTrajectory(gP, _keypoint2);

        //controlling motion along time
        Vector3 position = firstTrajectory * (1 - gP) + secondTrajectory * gP;

        return position;
    }

    public Quaternion finalOrientation(float time, Vector3 position)
    {
        float gF = findXOverTime(_time0, _time1, time);

        Quaternion qT1 = _new Toricmanifold(finalPosition.x,finalPosition.y, finalPosition.z, _keypoint1._target1, _keypoint1._target2).ComputeOrientation();
        Quaternion qT2 = _new Toricmanifold(finalPosition.x, finalPosition.y, finalPosition.z, _keypoint2._target1, _keypoint2._target2).ComputeOrientation();


        //controlling motion along time
        Quaternion rotation = qT1 * (1 - gF) + qT2 * gF;

        return rotation;
    }

    private Vector3 calculateTrajectory(float x, Toricmanifold targetAB)
    {
        Vector3 p0 = ToricComputing.FromWorldPosition(_keypoint1.toWorldPosition(), targetAB._target1, targetAB._target2);
        Vector3 p1 = ToricComputing.FromWorldPosition(_keypoint2.toWorldPosition(), targetAB._target1, targetAB._target2);

        Vector3[] targets = { targetAB._target1, targetAB._target2 };
        Vector3 AB =  targets[1] - targets[0];

        float alphaInterpolated = p0.x * x + p1.x * (1 - x);

        Vector3 vantageA0 = _keypoint1.toWorldPosition() - _keypoint1._target1;
        float distanceA0 = vantageA0.magnitude;
        vantageA0 = vantageA0.normalized;

        Vector3 vantageA1 = _keypoint2.toWorldPosition() - _keypoint1._target1;
        float distanceA1 = vantageA1.magnitude;
        vantageA0 = vantageA1.normalized;
        Vector3 interpolatedVa = x * vantageA0 + (1 - x) * vantageA1;
        float interpolatedDistanceA = x * distanceA0 + (1 - x) * distanceA1;

        Vector3 vantageB0 = _keypoint1.toWorldPosition() - _keypoint1._target2;
        float distanceB0 = vantageB0.magnitude;
        vantageB0 = vantageB0.normalized;

        Vector3 vantageB1 = _keypoint2.toWorldPosition() - _keypoint1._target2;
        float distanceB1 = vantageB1.magnitude;
        vantageB1 = vantageB1.normalized;
        Vector3 interpolatedVb = x * vantageB0 + (1 - x) * vantageB1;
        float interpolatedDistanceB = x * distanceB0 + (1 - x) * distanceB1;

        float thetaInterpolatedA = 2 * Vector3.Angle(interpolatedVa, AB);
        float thetaInterpolatedB = 2 * (Mathf.PI - Vector3.Angle(interpolatedBa, AB) - alphaInterpolated);


        Vector3 n = Vector3.Cross(interpolatedVa, AB).normalized;
        Vector3 z;

        Vector2 n2 = new Vector2(AB.x, AB.z);

        float tmp = n2[0];
        n2[0] = n2[1];
        n2[1] = -tmp;

        z = new Vector3(n2.x, 0, n2.y).normalized;
        float phiInterpolatedA = Vector3.SignedAngle(n, -AB, z) * Mathf.Deg2Rad - Mathf.PI / 2;

        Toricmanifold intersectionTa = new Toricmanifold( alphaInterpolated, thetaInterpolatedA, phiInterpolatedA, targets[0], targets[1] );
        float distanceAT = intersectionTa.toWorldPosition() - targets[0]; 


        Vector3 n = Vector3.Cross(interpolatedVb, AB).normalized;
        Vector3 z;

        Vector2 n2 = new Vector2(AB.x, AB.z);

        float tmp = n2[0];
        n2[0] = n2[1];
        n2[1] = -tmp;

        z = new Vector3(n2.x, 0, n2.y).normalized;
        float phiInterpolatedB = Vector3.SignedAngle(n, -AB, z) * Mathf.Deg2Rad - Mathf.PI / 2;

        Toricmanifold intersectionTb = new Toricmanifold(alphaInterpolated, thetaInterpolatedB, phiInterpolatedB, targets[0], targets[1]);
        float distanceBT = intersectionTb.toWorldPosition() - targets[1];


        float lambdaA = Mathf.Sin(Vector3.Angle(AB, interpolatedVa));
        float lambdaB = Mathf.Sin(Vector3.Angle(AB, interpolatedVb));

        return 0.5 * (targets[0] + targets[1] + lambdaA * ((interpolatedDistanceA + distanceAT * lambdaA) / (1 + lambdaA)) + lambdaB * ((interpolatedDistanceB + distanceBT * lambdaB) / (1 + lambdaB)));
    }

    private float findXOverTime(float time0, float time1, float time)
    {
        return (time - time0) / time1;
    }
}

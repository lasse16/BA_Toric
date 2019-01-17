﻿using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/**
* Basic test script for two targets
* 
*/
public class testScript : MonoBehaviour
{
    [Range(1, 179)]
    public float alpha;
    private float priorAlpha;

    [Range(0, 360)]
    public float theta;
    private float priorTheta;

    [Range(-180, 180)]
    public float phi;
    private float priorPhi;

    [Range(-90, 90)]
    public float tilt;
    private float priorTilt;

    public GameObject target1;
    public GameObject target2;

    [Tooltip("a Vector with the [min distance; max Distance] to Target A")]
    public Vector2 distanceToA;

    [Tooltip("a Vector with the [min distance; max Distance] to Target B")]
    public Vector2 distanceToB;

    [Tooltip("Desired screen position for Target A")]
    public Vector2 screenPos1;
    [Tooltip("Desired screen position for Target B")]
    public Vector2 screenPos2;

    public Vector2 sizeToReachA;
    public Vector2 sizeToReachB;

    public Vector3 vantageDirectionA;
    public float deviationAngleA;
    private Vector3 priorVantageA;

    public Vector3 vantageDirectionB;
    public float deviationAngleB;
    private Vector3 priorVantageB;

    public Vector2 visibilityInterval;
    public float samplingRate;

    public int time1;
    public int time2;
    /**
     * Creates either a cube or a camera at the calculated position with the given location
     * 
     */
    public void Start()
    {
        //StartDebug();
        //StartCamera();


        //testProjectedSize(); //TODO
        //testGetAlphaFromDistance();
        //testIntervalFromOnscreenPos(); //TODO
        testVantageAngleConstraint();
        //testVantageAngleConstraintA();
        //testVantageAngleConstraintB();
        //testFromWorld();
         //testInterpolation(); //TODO
        //testVisibility();
       //testAllConstraints();
        //visualizeTheta();
        //visualizePhi();
        //visualizeToricSpace();
        priorAlpha = alpha;
        priorTheta = theta;
        priorPhi = phi;
        priorTilt = tilt;
        priorVantageA = vantageDirectionA;
        priorVantageB = vantageDirectionB;







    }

    private void testAllConstraints()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Toricmanifold tm = tc.FinalConstraintCombination(samplingRate, vantageDirectionA, deviationAngleA, vantageDirectionB, deviationAngleB, distanceToA, distanceToB, screenPos1, screenPos2, visibilityInterval);

        //DEBUG
        Debug.Log("possiblePosition: " + tm.ToString());
        alpha = tm.getAlpha();
        theta = tm.getTheta();
        phi = tm.getPhi();

        StartCamera();
    }

    private void testVantageAngleConstraint()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Dictionary<float, Interval> phisA = tc.getThetaIntervallFromVantageBothTargets(vantageDirectionA, deviationAngleA, vantageDirectionB, deviationAngleB, 1 / samplingRate, true);


        //DEBUG
        float[] keys = new float[phisA.Keys.Count];
        phisA.Keys.CopyTo(keys, 0);

        phi = keys[UnityEngine.Random.Range(0, keys.Length - 2)];
        Interval betaRange;
        phisA.TryGetValue(phi, out betaRange);
        phi = phi * Mathf.Rad2Deg;
        Debug.Log(betaRange);
        theta = betaRange.getRandom() * Mathf.Rad2Deg;
        Debug.Log(theta);

        Debug.Log(phi);


        StartDebug();
    }

    private void testVantageAngleConstraintA()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Dictionary<float, Interval> phisA = tc._vantageCons.getPositionFromVantageOneTarget(target1.transform.position, vantageDirectionA, deviationAngleA);


        //DEBUG
        float[] keys = new float[phisA.Keys.Count];
        phisA.Keys.CopyTo(keys, 0);

        phi = keys[UnityEngine.Random.Range(0, keys.Length - 2)];
        Interval betaRange;
        phisA.TryGetValue(phi, out betaRange);
        phi = phi * Mathf.Rad2Deg;
        Debug.Log(betaRange);
        theta = 2 * betaRange.getRandom() * Mathf.Rad2Deg;
        Debug.Log(theta);

        Debug.Log(phi);


        StartDebug();
    }

    
    private void testInterpolation()
    {
        Toricmanifold test = new Toricmanifold(alpha, theta, phi, target1, target2);
        test.SetDesiredPosition(screenPos1, screenPos2);
        test.visualize(Color.green);
        Toricmanifold test2 = new Toricmanifold(alpha , theta, phi + 30, target1, target2);
        test2.SetDesiredPosition(screenPos1, screenPos2);
        test2.visualize(Color.red);

        ViewpointInterpolation interpol = new ViewpointInterpolation(test, test2, time1, time2);

        for (int i = time1; i < time2; i++)
        {
            Vector3 position = interpol.finalPosition(i);
            Quaternion rotation = interpol.finalOrientation(i, position);

            GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);
    
            cube.transform.position = position;
            cube.transform.rotation = rotation;
            cube.name = i + "";
        }
    }

    private void testFromWorld()
    {
        Toricmanifold test = new Toricmanifold(alpha, theta, phi, target1, target2);
        test.SetDesiredPosition(screenPos1, screenPos2);
        test.visualize(Color.red);
        Vector3 worldPos = test.ToWorldPosition();

        Vector3 toricRep = ToricComputing.FromWorldPosition(worldPos, target1.transform.position, target2.transform.position);
        Toricmanifold test2 = new Toricmanifold(toricRep.x * Mathf.Rad2Deg, toricRep.y * Mathf.Rad2Deg, toricRep.z * Mathf.Rad2Deg, target1, target2);
        test2.SetDesiredPosition(screenPos1, screenPos2);
        test2.visualize(Color.green); 
    }


    private void testVantageAngleConstraintB()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Dictionary<float, Interval> phisA = tc._vantageCons.getPositionFromVantageOneTarget(target2.transform.position, vantageDirectionB, deviationAngleB);


        //DEBUG
        float[] keys = new float[phisA.Keys.Count];
        phisA.Keys.CopyTo(keys, 0);

        phi = keys[UnityEngine.Random.Range(0, keys.Length - 2)];
        Interval betaRange;
        phisA.TryGetValue(phi, out betaRange);
        phi = phi * Mathf.Rad2Deg;
        Debug.Log(betaRange);
        theta = 2 * (Mathf.PI - alpha * Mathf.Deg2Rad - betaRange.getRandom()) * Mathf.Rad2Deg;
        Debug.Log(theta);

        Debug.Log(phi);


        StartDebug();
    }

    private void testProjectedSize()
    {


        ToricComputing tc = new ToricComputing(target1, target2);
        Vector2 projectedSizeA = tc.DistanceFromProjectedSize(sizeToReachA, 0.5f, 1,true).ToVector();
        Vector2 projectedSizeB = tc.DistanceFromProjectedSize(sizeToReachB, 0.5f, 2,true).ToVector();
        //DEBUG
        Debug.Log(projectedSizeA + ";" + projectedSizeB);
        Dictionary<float, Interval> alphas = tc._alphaComputer.getIntervalOfAcceptedAlpha(projectedSizeA, projectedSizeB);

        foreach (KeyValuePair<float, Interval> a in alphas)
        {
            Debug.Log(a.ToString());
        }

        //DEBUG
        float[] keys = new float[alphas.Keys.Count];
        alphas.Keys.CopyTo(keys, 0);
        theta = keys[UnityEngine.Random.Range(0, alphas.Keys.Count - 2)];
        Interval alphaRange;
        alphas.TryGetValue(theta, out alphaRange);
        alpha = alphaRange.getRandom();
        Debug.Log(theta);
        Debug.Log(alphaRange);
        Debug.Log(alpha);



        StartDebug();
    }

    private void StartCamera()
    {
        Toricmanifold test = new Toricmanifold(alpha, theta, phi, target1, target2);
        test.SetDesiredPosition(screenPos1, screenPos2);

        Camera _main = Camera.main;

        Vector3 posTest = test.ToWorldPosition();
        Quaternion rotTest = test.ComputeOrientation( tilt);

        _main.transform.position = posTest;
        _main.transform.rotation = rotTest;
    }

    private void StartDebug()
    {
        Toricmanifold test = new Toricmanifold(alpha, theta, phi, target1, target2);
        test.SetDesiredPosition(screenPos1, screenPos2);


        Vector3 posTest = test.ToWorldPosition();
        Quaternion rotTest = test.ComputeOrientation( tilt);
        GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);

        cube.transform.position = posTest;
        cube.transform.rotation = rotTest;
        cube.name = alpha + " : " + theta + " : " + phi;


    }

    private void testGetAlphaFromDistance()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Dictionary<float, Interval> alphaInv = tc._alphaComputer.getIntervalOfAcceptedAlpha(distanceToA, distanceToB,0.05f,true);

        Interval alphaRange;
        float[] keys = new float[alphaInv.Keys.Count];
        alphaInv.Keys.CopyTo(keys, 0);
        if (alphaInv.Keys.Count - 2 <= 0) throw new Exception("No intersection");
        theta = keys[UnityEngine.Random.Range(0, alphaInv.Keys.Count - 2)];
        alphaInv.TryGetValue(theta, out alphaRange);


        alpha = alphaRange.getRandom() * Mathf.Rad2Deg;
        theta *= Mathf.Rad2Deg;

        Debug.Log(alphaRange);
        Toricmanifold test = new Toricmanifold(alpha, theta, phi, target1, target2);
        Vector3 posTest = test.ToWorldPosition();
        
        GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);

        //DEBUG
        Debug.Log("distance: " + (posTest - target1.transform.position).magnitude);
        Debug.Log("distance: " + (posTest - target2.transform.position).magnitude);

        cube.transform.position = posTest;
        
    }

    

    private void testIntervalFromOnscreenPos()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Interval alphaTest = tc._alphaComputer.getAlphaIntervalFromOnscreenPositions(screenPos1, screenPos2,true);



        alpha = alphaTest.getRandom();

        StartCamera();

    }


    

    private void visualizeTheta()
    {
        float givenTheta = theta;
        Toricmanifold tm = new Toricmanifold(alpha, 1f, phi, target1, target2);
        new Ellipse(.01f, .01f, tm.ToWorldPosition(), Vector3.up).Draw(Color.red);

        float maxTheta = tm.getMaxTheta();

        List<Vector3> thetaPos = new List<Vector3>();

        for (int i = 2; i < maxTheta; i++)
        {
            Toricmanifold manifold = new Toricmanifold(alpha, i, phi, target1, target2);
            thetaPos.Add(manifold.ToWorldPosition());
        }

        Vector3 priorPos = tm.ToWorldPosition();
        foreach (Vector3 pos in thetaPos)
        {
            Debug.DrawLine(priorPos, pos, Color.gray, Mathf.Infinity, false);
            priorPos = pos;
        }

        tm = new Toricmanifold(alpha, givenTheta, phi, target1, target2);
        new Ellipse(.01f, .01f, tm.ToWorldPosition(), Vector3.up).Draw(Color.blue);
    }

    private void visualizePhi()
    {
        Toricmanifold tm = new Toricmanifold(alpha, theta, 1, target1, target2);
        int maxPhi = 180;

        List<Vector3> thetaPos = new List<Vector3>();

        for (int i = -maxPhi; i < maxPhi; i++)
        {
            Toricmanifold manifold = new Toricmanifold(alpha, theta, i, target1, target2);
            thetaPos.Add(manifold.ToWorldPosition());
        }

        Vector3 priorPos = tm.ToWorldPosition();
        foreach (Vector3 pos in thetaPos)
        {
            Debug.DrawLine(priorPos, pos, Color.grey, Mathf.Infinity, false);
            priorPos = pos;
        }
    }

    private void visualizeToricSpace()
    {

        Toricmanifold tm = new Toricmanifold(alpha, 1, 1, target1, target2);
        int maxPhi = 180;
        float maxTheta = tm.getMaxTheta();


        for (int j = 0; j < maxTheta; j += 3)
        {
            List<Vector3> phiPos = new List<Vector3>();

            for (int i = -maxPhi; i < maxPhi; i++)
            {
                Toricmanifold manifold = new Toricmanifold(alpha, j, i, target1, target2);
                phiPos.Add(manifold.ToWorldPosition());
            }

            Vector3 priorPos = tm.ToWorldPosition();
            foreach (Vector3 pos in phiPos)
            {
                Debug.DrawLine(priorPos, pos, Color.grey, Mathf.Infinity, false);
                priorPos = pos;
            }
        }
    }


    void Update()
    {
        if (priorAlpha != 0 && priorTheta != 0)
        {
            if (priorAlpha != alpha || priorPhi != phi || priorTheta != theta || priorTilt != tilt || !priorVantageA.Equals(vantageDirectionA) || !priorVantageB.Equals(vantageDirectionB)) Start();
        }
    }
}
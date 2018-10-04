using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/**
* Basic test script for two targets
* 
*/
public class testScript : MonoBehaviour
{
    [Range(1,179)]
    public float alpha;
    private float priorAlpha;

    [Range (0,360)]
    public float theta;
    private float priorTheta;

    [Range(-180, 180)]
    public float phi;
    private float priorPhi;

    [Range(-90,90)]
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

    /**
     * Creates either a cube or a camera at the calculated position with the given location
     * 
     */
    public void Start()
    {
        //StartDebug();
        //StartCamera();
        //StartComputing();
        // testGetAlphaFromDistanceB(); //TODO
        //testIntervalFromOnscreenPos();
        //testIntervalFromB();
       testVantageAngleConstraintA(); //TODO
        //testVisibility();
       // testAllConstraints(); //TODO
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
        Toricmanifold tm = tc.FinalConstraintCombination( samplingRate, vantageDirectionA,  deviationAngleA,  vantageDirectionB,  deviationAngleB,  distanceToA,  distanceToB,  screenPos1,  screenPos2,  visibilityInterval);

        //DEBUG
        Debug.Log("possiblePosition: " + tm.ToString());
        alpha = tm.getAlpha();
        theta = tm.getTheta();
        phi = tm.getPhi();

        StartCamera();
    }

    private void testVantageAngleConstraintA()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Dictionary<float, Interval> phisA = tc.getPositionFromVantageOneTarget(1, vantageDirectionA, deviationAngleA);
        

        //DEBUG
        float[] keys = new float[phisA.Keys.Count];
        phisA.Keys.CopyTo(keys, 0);
        
        phi = keys[UnityEngine.Random.Range(0, keys.Length - 2)];
        Interval thetaRange;
        phisA.TryGetValue(phi, out thetaRange);
        phi = phi * Mathf.Rad2Deg;
        Debug.Log(thetaRange);
        theta = thetaRange.getRandom() * Mathf.Rad2Deg;
        Debug.Log(theta);
        
        Debug.Log(phi);


        StartDebug();
    }

    private void StartComputing()
    {


        ToricComputing tc = new ToricComputing(target1, target2);
        Vector2 projectedSizeA = tc.DistanceFromProjectedSize(sizeToReachA, 0.5f, target1).toVector();
        Vector2 projectedSizeB = tc.DistanceFromProjectedSize(sizeToReachB, 0.5f, target1).toVector();
        //DEBUG
        Debug.Log(projectedSizeA + ";" + projectedSizeB);
        Dictionary<float, Interval> alphas = tc.getIntervalOfAcceptedAlpha(projectedSizeA, projectedSizeB);

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
        Quaternion rotTest = test.ComputeOrientation(posTest, tilt);

        _main.transform.position = posTest;
        _main.transform.rotation = rotTest;
    }

    private void StartDebug()
    {
        Toricmanifold test = new Toricmanifold(alpha, theta, phi, target1, target2);
        test.SetDesiredPosition(screenPos1, screenPos2);

        Vector3 posTest = test.ToWorldPosition();
        Quaternion rotTest = test.ComputeOrientation(posTest, tilt);
        GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);

        cube.transform.position = posTest;
        cube.transform.rotation = rotTest;
        cube.name = alpha + " : " + theta + " : " + phi;

        
    }

    private void testGetAlphaFromDistanceA()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        float alphaTestMin = tc.testDistanceFromA(distanceToA[0], theta);
        float alphaTestMax = tc.testDistanceFromA(distanceToA[1], theta);

        Toricmanifold test = new Toricmanifold(alphaTestMax, theta, phi, target1, target2);
        Vector3 posTest = test.ToWorldPosition();
        Quaternion rotTest = test.ComputeOrientation(posTest, tilt);
        GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);

        //DEBUG
        Debug.Log("distance: " + (posTest - target1.transform.position).magnitude);

        cube.transform.position = posTest;
        cube.transform.rotation = rotTest;
    }

    private void testGetAlphaFromDistanceB()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Interval alphaTestMin = tc.GetAlphaFromDistanceToB(distanceToB[0], theta);
        Interval alphaTestMax = tc.GetAlphaFromDistanceToB(distanceToB[1], theta);


        Interval alphaInterval =  alphaTestMin - alphaTestMax ;
        alpha = alphaInterval.getRandom();

        //DEBUG
        Debug.Log(alphaInterval);


        StartDebug();
    }

    private void testIntervalFromOnscreenPos()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Interval alphaTest = tc.getAlphaIntervalFromOnscreenPositions(screenPos1, screenPos2);

        

        alpha = alphaTest.getRandom();

        StartCamera();

    }

    private void testIntervalFromB()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Dictionary <float,Interval> alphas = tc.getIntervalFromB(screenPos2.x,screenPos2.y);

       
        float[] keys = new float[alphas.Keys.Count];
        alphas.Keys.CopyTo(keys, 0);
        theta = keys[UnityEngine.Random.Range(0, alphas.Keys.Count - 2)];
        Interval alphaRange;
        alphas.TryGetValue(theta, out alphaRange);
        alpha = alphaRange.getRandom();


        StartCamera();
    }

    private void testVisibility()
    {
        ToricComputing tc = new ToricComputing(target1, target2);
        Toricmanifold tm = new Toricmanifold(alpha, theta, phi, target1, target2);
        Debug.Log(tc.testVisibility(tm));

        StartCamera();
    }


    void Update()
    {
        if (priorAlpha != alpha || priorPhi != phi || priorTheta != theta || priorTilt != tilt || !priorVantageA.Equals(vantageDirectionA) || !priorVantageB.Equals(vantageDirectionB)) Start();
    }
}
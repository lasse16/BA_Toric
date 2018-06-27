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
    public float alpha;
    public float theta;
    public float phi;
    public float tilt;
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

    /**
     * Creates either a cube or a camera at the calculated position with the given location
     * 
     */
    public void Start()
    {
        //StartDebug();
        //StartCamera();
        StartDistance();
       

       

       
        

    }

    private void StartDistance()
    {
        //CHECK FOR EXTRA CLASS

        ToricComputing tc = new ToricComputing(target1, target2);
        Dictionary<float, Intervall> alphas =  tc.getIntervalOfAcceptedAlpha(distanceToA,distanceToB);

       

        Intervall outInt;
        float[] thetas = new float[alphas.Keys.Count];
        alphas.Keys.CopyTo(thetas, 0);
        theta = thetas[UnityEngine.Random.Range(0, thetas.Length)];

        //DEBUG
        foreach (float t in thetas)
        {
            Intervall output;
            alphas.TryGetValue(t, out output);
            Debug.Log(t + ":" + output.Length() + "-");
        }       
        //DEBUG
        Debug.Log(theta);
        Debug.Log(ToricComputing.FloatArrayToString(thetas));

        alphas.TryGetValue(theta, out outInt);
        alpha = outInt.getRandom();

        //DEBUG
        Debug.Log(outInt);
        Debug.Log(alpha);
        

        Toricmanifold test = new Toricmanifold(alpha, theta, phi, target1, target2);
        test.SetDesiredPosition(screenPos1, screenPos2);

        Vector3 posTest = test.ToWorldPosition();
        Quaternion rotTest = test.ComputeOrientation(posTest, tilt);
        GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);

        cube.transform.position = posTest;
        cube.transform.rotation = rotTest;
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
    }
}

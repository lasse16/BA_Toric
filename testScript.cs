using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class testScript : MonoBehaviour
{
    public float alpha;
    public float theta;
    public float phi;
    public float tilt;
    public GameObject target1;
    public GameObject target2;

    [Tooltip("Desired screen position for Target A")]
    public Vector2 screenPos1;
    [Tooltip("Desired screen position for Target B")]
    public Vector2 screenPos2;


    public void Start()
    {


        Toricmanifold test = new Toricmanifold(alpha, theta, phi, target1, target2);
        test.setDesiredPosition(screenPos1, screenPos2);

        Vector3 posTest = test.ToWorldPosition();
        Quaternion rotTest = test.ComputeOrientation(posTest, tilt);
        GameObject cube = GameObject.CreatePrimitive(PrimitiveType.Cube);

        Camera _main = Camera.main;

        //_main.transform.position = posTest;
        //_main.transform.rotation = rotTest;
        cube.transform.position = posTest;
        cube.transform.rotation = rotTest;

    }
}

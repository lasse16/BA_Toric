using System.Collections.Generic;
using UnityEngine;

internal class Cone
{
    private Vector3 _apex;
    private Vector3 _prefferedVantageAngle;
    //allowed deviation around the preffered vantage angle (degrees) 
    private float _deviationAngle;


    public Cone(Vector3 prefferedVantageAngle, float deviationAngle, Vector3 apex)
    {
        _prefferedVantageAngle = prefferedVantageAngle;
        _deviationAngle = deviationAngle;
        _apex = apex;
    }

    public Cone(Vector3 prefferedVantageAngle, float deviationAngle): this(prefferedVantageAngle,deviationAngle,Vector3.zero)
    {    }

    public Vector3 getMidLineOfCone()
    {
        return _prefferedVantageAngle;
    }

    public float getDiviationAngle()
    {
        return _deviationAngle;
    }

    public bool IsInside(Vector3 testVector)
    {
        return Vector3.Angle(testVector, _prefferedVantageAngle) <= _deviationAngle;
    }


    public void draw(float coneHeight)
    {
       
            Vector3[] positions;
            int resolution = 16;
            float a = Mathf.Tan(_deviationAngle) * coneHeight;
           



            positions = new Vector3[resolution + 1];
            float b = Vector3.Angle(Vector3.up, _prefferedVantageAngle)- 90;
            Quaternion q = Quaternion.LookRotation( _prefferedVantageAngle);
            Vector3 center = _apex + _prefferedVantageAngle * coneHeight;

            for (int i = 0; i <= resolution; i++)
            {
                float angle = (float)i / (float)resolution * 2.0f * Mathf.PI;
                positions[i] = new Vector3(a * Mathf.Cos(angle), a * Mathf.Sin(angle), 0.0f);
                positions[i] = q * positions[i] + center;
            }

            Vector3 priorPoint = center;
            foreach (Vector3 point in positions)
            {
                Debug.DrawLine(priorPoint, point, Color.green, Mathf.Infinity, false);
                Debug.DrawLine(_apex, point, Color.green, Mathf.Infinity, false);
            priorPoint = point;
            }

        

    }




}